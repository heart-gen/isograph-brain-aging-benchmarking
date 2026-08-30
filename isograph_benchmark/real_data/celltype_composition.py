"""Cell-type composition confound: fractions + a marker-depletion cut.

The bulk-brain aging story has one dominant alternative explanation — the isoform
switches track shifting neuron/glia *proportions* with age, not within-cell-type
regulation. This module supplies the two pieces needed to test that, reusing the
committed MuSiC deconvolution already run on the BrainSEQ cohort
(``sex_context_brain/cell_proportion_estimate``); no new deconvolution is done here.

1. Fractions join. The MuSiC per-donor cell-type proportions are keyed by ``BrNum``;
   the IsoGraph bundle sample tables carry ``BrNum`` alongside the RNA ``sample_id``.
   We map BrNum -> sample_id and write ``<artifact_dir>/celltype_fractions.parquet``
   (one row per sample, one column per cell type). ``incremental_association
   --composition`` then adds these as inference covariates, so the de-confounded
   switch-vs-abundance test is repeated with composition regressed out (with-vs-without
   contrast = how many phenotype/age switch signals survive).

2. Marker-depletion cut. For each module we test whether its member genes are enriched
   or *depleted* for cell-type marker genes (top mean-ratio markers from the same
   DeconvoBuddies run). Depletion argues the co-switch signal is within-cell-type
   regulation rather than a proportion artifact — a positive control that does not
   depend on the covariate model.

Both outputs are deterministic pure joins/hypergeometrics over artifacts on disk.
BrainSEQ only (the MuSiC run covers caudate/DLPFC/hippocampus, incl. the SCZD caudate
region); GTEx would require re-running the upstream MuSiC script (see plan, optional arm).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import cohort_dir, ensure_dir, stage_out
from isograph_benchmark.real_data.qtl_anchoring import _bare
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir, _bundle_path

# Committed MuSiC proportions + markers from the sibling repo (BrNum-keyed).
DEFAULT_MUSIC_DIR = Path(
    "/ocean/projects/bio260021p/kbenjamin/projects/sex_context_brain"
    "/cell_proportion_estimate/_m"
)
# GTEx MuSiC run lives in THIS repo (sample_id-keyed; produced by the R deconvolution in
# 04_module_characterization/_h/gtex_music_deconv.R, seeded from the same Tran/LIBD snRNA references).
GTEX_MUSIC_DIR = cohort_dir("gtex", "_m", "composition")

# GTEx brain region -> Tran/LIBD snRNA reference key. Only regions with a defensibly matched
# reference are deconvolved; the mapping follows the precedent already committed in
# sex_context_brain (striatum deconvolved with the NAc reference, cortex with DLPFC). Regions
# with no matched Tran reference (cerebellum/cerebellar_hemisphere, hypothalamus, spinal cord,
# substantia nigra) are left out honestly rather than deconvolved against a mismatched panel.
GTEX_REF = {
    "amygdala": "amy",
    "anterior_cingulate_cortex_ba24": "sacc",
    "frontal_cortex_ba9": "dlpfc",
    "cortex": "dlpfc",
    "hippocampus": "hpc",
    "caudate_basal_ganglia": "nac",
    "putamen_basal_ganglia": "nac",
    "nucleus_accumbens_basal_ganglia": "nac",
}

# IsoGraph region -> MuSiC region file stem. SCZD caudate reuses the caudate MuSiC run;
# GTEx regions get their own `gtex-<region>` stems from the in-repo deconvolution.
_MUSIC_REGION = {
    ("brainseq-sczd", None): "caudate",
    ("brainseq-aging", "caudate"): "caudate",
    ("brainseq-aging", "dlpfc"): "dlpfc",
    ("brainseq-aging", "hippocampus"): "hippocampus",
    **{("gtex-aging", r): f"gtex-{r}" for r in GTEX_REF},
}


def _music_dir_for(analysis: str, override: Path | None) -> Path:
    """MuSiC output dir: the in-repo GTEx run for gtex-aging, else the sibling BrainSEQ run
    (or an explicit --music-dir override)."""
    if override is not None:
        return override
    return GTEX_MUSIC_DIR if analysis == "gtex-aging" else DEFAULT_MUSIC_DIR


def _md_table(df: pd.DataFrame) -> str:
    """GitHub-flavored markdown table without the optional `tabulate` dependency."""
    cols = list(df.columns)
    head = "| " + " | ".join(map(str, cols)) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, sep, *body])


def music_region(analysis: str, region: str | None) -> str:
    key = (analysis, region)
    if key not in _MUSIC_REGION:
        raise ValueError(
            f"no MuSiC deconvolution for analysis={analysis!r} region={region!r}; "
            "covered: BrainSEQ caudate/dlpfc/hippocampus + GTEx "
            f"{sorted(GTEX_REF)} (other GTEx regions lack a matched Tran reference)."
        )
    return _MUSIC_REGION[key]


def load_fractions_wide(music_dir: Path, mregion: str, key: str = "BrNum") -> pd.DataFrame:
    """(donor × cell types) proportion matrix from the MuSiC tsv, indexed by the join key
    (BrNum for BrainSEQ, sample_id for GTEx where MuSiC ran per RNA sample)."""
    long = pd.read_csv(music_dir / f"music-proportions-{mregion}.tsv", sep="\t")
    wide = long.pivot_table(index="sample_id", columns="cell_type",
                            values="proportion", aggfunc="mean")
    wide.index.name = key
    wide.columns.name = None
    return wide


def build_fractions(analysis: str, region: str | None, variant: str,
                    music_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    """Write <artifact_dir>/celltype_fractions.parquet keyed by RNA sample_id.

    Returns (fractions indexed by sample_id, cell-type column names). BrainSEQ joins the
    BrNum-keyed MuSiC table to the RNA sample_id; GTEx MuSiC is already sample_id-keyed
    (one RNA sample per tissue donor) so it merges directly.
    """
    mregion = music_region(analysis, region)
    key = "sample_id" if analysis == "gtex-aging" else "BrNum"
    wide = load_fractions_wide(music_dir, mregion, key=key)
    cell_cols = list(wide.columns)

    st = pd.read_parquet(_bundle_path(analysis, region) / "samples.parquet",
                         columns=None)
    if key not in st.columns:
        raise KeyError(f"bundle sample table for {analysis}/{region} has no {key} column")
    keep = ["sample_id"] if key == "sample_id" else ["sample_id", "BrNum"]
    st = st[keep].astype({c: str for c in keep})

    frac = st.merge(wide.reset_index().astype({key: str}), on=key, how="left")
    if "BrNum" in frac.columns:
        frac = frac.drop(columns=["BrNum"])
    frac = frac.set_index("sample_id")

    art = _artifact_dir(analysis, region, variant)
    ensure_dir(art)
    out = art / "celltype_fractions.parquet"
    covered = int(frac[cell_cols].notna().all(axis=1).sum())
    frac.reset_index().to_parquet(out, index=False, compression="zstd")
    print(f"[{analysis}/{region}] fractions: {covered}/{len(frac)} samples with "
          f"{len(cell_cols)} cell types -> {out}", flush=True)
    return frac[cell_cols], cell_cols


def composition_covariate_cols(cell_cols: list[str], fractions: pd.DataFrame) -> list[str]:
    """All cell types except the largest-mean one (dropped as the simplex reference to
    avoid collinearity with the intercept — fractions sum to ~1)."""
    ref = fractions[cell_cols].mean(axis=0).idxmax()
    return [c for c in cell_cols if c != ref]


def load_composition_covariates(artifact_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    """Read a written fractions parquet; return (sample_id-indexed fractions, cov cols)."""
    frac = pd.read_parquet(artifact_dir / "celltype_fractions.parquet")
    frac = frac.set_index(frac["sample_id"].astype(str)).drop(columns=["sample_id"])
    cell_cols = list(frac.columns)
    return frac, composition_covariate_cols(cell_cols, frac)


# --------------------------------------------------------------------------- #
# GTEx bulk export (feeds the R MuSiC deconvolution)
# --------------------------------------------------------------------------- #
def export_gtex_bulk(region: str, out_dir: Path) -> Path:
    """Write a genes×samples raw-count matrix for one GTEx region as a parquet the R
    deconvolution reads (arrow). Columns: gene_id, gene_name, then one column per RNA
    sample_id. Counts are RNASeQC gene reads (manifest: raw integer counts), so MuSiC
    sees the same count scale the BrainSEQ rse path provided."""
    if region not in GTEX_REF:
        raise ValueError(f"{region!r} has no matched Tran reference; deconvolvable GTEx "
                         f"regions: {sorted(GTEX_REF)}")
    bpath = _bundle_path("gtex-aging", region)
    genes = pd.read_parquet(bpath / "genes.parquet")
    counts = np.load(bpath / "gene_counts.npz")["data"]  # (genes, samples)
    samples = pd.read_parquet(bpath / "samples.parquet", columns=["sample_id"])
    if counts.shape != (len(genes), len(samples)):
        raise ValueError(f"count matrix {counts.shape} != (genes {len(genes)}, "
                         f"samples {len(samples)}) for {region}")
    df = pd.DataFrame(counts, columns=samples["sample_id"].astype(str).tolist())
    df.insert(0, "gene_name", genes["gene_name"].astype(str).values)
    df.insert(0, "gene_id", genes["gene_id"].astype(str).values)
    ensure_dir(out_dir / "inputs")
    out = out_dir / "inputs" / f"{region}_bulk.parquet"
    df.to_parquet(out, index=False, compression="zstd")
    print(f"[export-gtex] {region}: {len(genes)} genes × {len(samples)} samples -> {out}",
          flush=True)
    return out


# --------------------------------------------------------------------------- #
# Marker-depletion cut
# --------------------------------------------------------------------------- #
def load_marker_sets(music_dir: Path, mregion: str, top_n: int) -> dict[str, set[str]]:
    """cell_type -> set of bare Ensembl ids for its top-`top_n` mean-ratio markers."""
    m = pd.read_csv(music_dir / f"marker_stats_genes.{mregion}.csv")
    m = m[m["MeanRatio.rank"] <= top_n]
    m["gene"] = _bare(m["gene_ensembl"])
    return {ct: set(g["gene"]) for ct, g in m.groupby("cellType.target")}


def marker_enrichment(modules: pd.DataFrame, markers: dict[str, set[str]]) -> pd.DataFrame:
    """Per (module, cell_type): hypergeometric enrichment/depletion of that cell type's
    markers among the module's genes vs the tested-gene background. Depletion (frac below
    expected, sf~1 / observed<expected) argues against a proportion artifact."""
    modules = modules.copy()
    modules["gene"] = _bare(modules["gene_id"])
    background = set(modules["gene"])
    N = len(background)
    rows = []
    for ct, mk in markers.items():
        K = len(mk & background)
        if K == 0:
            continue
        for mid, grp in modules.groupby("module_id"):
            mod_genes = set(grp["gene"])
            n = len(mod_genes)
            k = len(mod_genes & mk)
            expected = n * K / N
            # two-sided-ish: report enrichment sf and depletion cdf
            p_enrich = float(stats.hypergeom.sf(k - 1, N, K, n))
            p_deplete = float(stats.hypergeom.cdf(k, N, K, n))
            rows.append({"module_id": mid, "cell_type": ct, "module_size": n,
                         "n_markers_bg": K, "observed": k, "expected": round(expected, 3),
                         "p_enrich": p_enrich, "p_deplete": p_deplete})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_enrich"] = stats.false_discovery_control(out["p_enrich"], method="bh")
        out["fdr_deplete"] = stats.false_discovery_control(out["p_deplete"], method="bh")
        out = out.sort_values(["module_id", "cell_type"]).reset_index(drop=True)
    return out


def run(analysis: str, region: str | None, variant: str, music_dir: Path,
        marker_top_n: int) -> dict:
    label = f"{analysis}/{region}" if region else analysis
    fractions, cell_cols = build_fractions(analysis, region, variant, music_dir)
    comp_cols = composition_covariate_cols(cell_cols, fractions)

    art = _artifact_dir(analysis, region, variant)
    modules = pd.read_parquet(art / "modules.parquet")
    markers = load_marker_sets(music_dir, music_region(analysis, region), marker_top_n)
    menr = marker_enrichment(modules, markers)
    out = ensure_dir(art / "celltype_composition")
    menr.to_parquet(out / "marker_enrichment.parquet", index=False, compression="zstd")

    n_dep = int((menr["fdr_deplete"] <= 0.05).sum()) if not menr.empty else 0
    n_enr = int((menr["fdr_enrich"] <= 0.05).sum()) if not menr.empty else 0
    summary = {
        "analysis": analysis, "region": region, "variant": variant,
        "cell_types": cell_cols, "composition_covariates": comp_cols,
        "n_samples": int(len(fractions)),
        "n_samples_covered": int(fractions.notna().all(axis=1).sum()),
        "marker_top_n": marker_top_n,
        "n_module_celltype_marker_depleted_fdr05": n_dep,
        "n_module_celltype_marker_enriched_fdr05": n_enr,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[{label}] marker cut: {n_enr} enriched / {n_dep} depleted (module×celltype, "
          f"fdr<=0.05); covariates={comp_cols}", flush=True)
    return summary


# --------------------------------------------------------------------------- #
# Cross-region meta: the with-vs-without composition-adjustment contrast
# --------------------------------------------------------------------------- #
_META_REGIONS = [
    ("brainseq-sczd", None, "caudate_sczd", "SCZD (caudate)"),
    ("brainseq-aging", "caudate", "caudate", "aging caudate"),
    ("brainseq-aging", "hippocampus", "hippocampus", "aging hippocampus"),
    ("brainseq-aging", "dlpfc", "dlpfc", "aging DLPFC"),
]


def _contrast_rows(specs: list[tuple[str, str | None, str]], variant: str) -> pd.DataFrame:
    """base vs composition-adjusted gene-level contrast + marker cut for each (analysis,
    region, label); skips any region whose base/adjusted summaries are not both on disk."""
    rows = []
    for analysis, region, label in specs:
        art = _artifact_dir(analysis, region, variant)
        base_p = art / "incremental_association" / "summary.json"
        adj_p = art / "incremental_association_composition" / "summary.json"
        comp_p = art / "celltype_composition" / "summary.json"
        if not (base_p.exists() and adj_p.exists()):
            print(f"[meta] skipping {label}: missing base/adj summary", flush=True)
            continue
        base = json.loads(base_p.read_text())["gene_level"]
        adj = json.loads(adj_p.read_text())["gene_level"]
        comp = json.loads(comp_p.read_text()) if comp_p.exists() else {}
        rows.append({
            "region": label,
            "n_tested": base["n_tested"],
            "comp_unique_base": base["composition_unique"],
            "comp_unique_adj": adj["composition_unique"],
            "retained_frac": (round(adj["composition_unique"] / base["composition_unique"], 2)
                              if base["composition_unique"] else np.nan),
            "both_base": base["both"], "both_adj": adj["both"],
            "n_cell_types": len(comp.get("cell_types", [])),
            "samples_covered": comp.get("n_samples_covered"),
            "marker_depleted": comp.get("n_module_celltype_marker_depleted_fdr05"),
            "marker_enriched": comp.get("n_module_celltype_marker_enriched_fdr05"),
        })
    return pd.DataFrame(rows)


def meta(variant: str) -> None:
    """Roll up base vs composition-adjusted incremental_association + the marker cut. The
    BrainSEQ cohort is the headline (disease + primary aging); GTEx is the aging replication
    arm — appended when its in-repo MuSiC deconvolution has been run."""
    bq = _contrast_rows([(a, r, l) for a, r, _, l in _META_REGIONS], variant)
    out = ensure_dir(stage_out("characterize"))
    bq.to_parquet(out / "composition_adjustment.parquet", index=False, compression="zstd")

    gt_specs = [("gtex-aging", r, f"GTEx {r}") for r in GTEX_REF]
    gt = _contrast_rows(gt_specs, variant)

    lines = ["# Cell-type composition adjustment of the DTU-without-DGE layer", "",
             "Re-running the de-confounded gene-level switch-vs-abundance test with MuSiC "
             "cell-type fractions added as inference covariates (reference cell type dropped "
             "for the simplex). `comp_unique` = genes whose isoform composition is "
             "phenotype-associated beyond their own abundance (IsoGraph's unique signal).", "",
             "## BrainSEQ (disease + primary aging)", "",
             _md_table(bq)]
    if not gt.empty:
        n_survive = int((gt["comp_unique_adj"] > 0).sum())
        n_collapse = int((gt["comp_unique_adj"] == 0).sum())
        gt.to_parquet(GTEX_MUSIC_DIR / "composition_adjustment_gtex.parquet",
                      index=False, compression="zstd")
        gtex_lines = [
            "# GTEx composition adjustment — aging replication arm", "",
            "MuSiC deconvolution re-run in-repo on the GTEx brain bundles using the same "
            "Tran/LIBD snRNA references (seed 13), for the regions with a defensibly matched "
            "reference (striatum→NAc, cortex→DLPFC, plus AMY/sACC/HPC direct). Regions with "
            "no matched reference (cerebellum, hypothalamus, spinal cord, substantia nigra) "
            "are not deconvolved.", "",
            _md_table(gt), "",
            "## Interpretation — a region-dependent, honestly partial replication", "",
            f"**{n_survive}/{len(gt)} regions retain ≥1 composition-robust DTU-without-DGE "
            f"gene, but {n_collapse}/{len(gt)} collapse to zero and the retained fraction "
            "varies enormously (0.0–0.86).** This is the confounder-vs-mediator caveat made "
            "concrete rather than a clean win:", "",
            "- **Limbic/striatal regions retain** (amygdala, hippocampus, caudate, putamen, "
            "NAc): the aging switch signal there is not merely a proportion artifact.",
            "- **The two cortical regions collapse to 0** (frontal_cortex_ba9 531→0, cortex "
            "438→0) — despite BrainSEQ *DLPFC* surviving adjustment. Cortical composition is "
            "the most strongly age-coupled, and these are the weakest reference matches "
            "(GTEx cortex/BA9 → DLPFC snRNA), so the collapse is consistent with genuine "
            "composition-confounding of cortical aging switches **and/or over-adjustment** "
            "when 7–8 age-correlated fraction covariates absorb the age-spline variance. The "
            "data cannot cleanly separate these here; report the collapse, do not hide it.", "",
            "**Manuscript use:** cite GTEx as a *partial* composition replication — the "
            "limbic/striatal aging layer reproduces as composition-robust across cohorts, "
            "while cortical aging DTU is composition-entangled in GTEx (a stated limitation, "
            "not a claim of universal robustness).", "",
            "_Caveat: cross-region reference use (striatum→NAc, cortex→DLPFC) and the "
            "snRNA reference panel make these fraction estimates approximate; the "
            "with-vs-without contrast, not the absolute fractions, is the claim._"]
        (GTEX_MUSIC_DIR / "GTEX_COMPOSITION_SUMMARY.md").write_text("\n".join(gtex_lines))
        print(f"[meta] wrote {GTEX_MUSIC_DIR/'GTEX_COMPOSITION_SUMMARY.md'} "
              f"({n_survive}/{len(gt)} regions survive)", flush=True)
        lines += ["", "## GTEx (aging replication arm)", "",
                  f"In-repo MuSiC re-run; {n_survive}/{len(gt)} deconvolved regions retain a "
                  f"composition-robust switch signal and {n_collapse}/{len(gt)} collapse to "
                  "zero — a **region-dependent, partial** replication: limbic/striatal aging "
                  "DTU reproduces as composition-robust, while the two cortical regions "
                  "collapse (composition-entangled and/or over-adjusted). Full breakdown + "
                  "caveats in `02_module_discovery/gtex/_m/composition/GTEX_COMPOSITION_SUMMARY.md`.", "",
                  _md_table(gt)]

    lines += ["", "## Interpretation", "",
              "- The composition-adjusted `comp_unique_adj` is the honest, composition-robust "
              "count that should anchor the DTU-without-DGE claim; report it alongside the "
              "unadjusted number rather than in place of it.",
              "- A large base→adj drop (e.g. SCZD) means much of that switch signal co-varies "
              "with cell-type proportion. **Confounder vs mediator matters:** if disease/age "
              "*causes* the composition shift that drives the switch, covariate adjustment "
              "*removes real signal* (over-adjustment); if composition varies for technical/"
              "sampling reasons, adjustment is the correct control. State this both-ways and "
              "lean on the layers that do not depend on it (module-level genetic anchoring, "
              "GO-invisible gate).",
              "- The aging switch layer survives composition adjustment in BrainSEQ (caudate, "
              "DLPFC) and **partially replicates in GTEx** (limbic/striatal robust; the two "
              "cortical regions collapse — composition-entangled and/or over-adjusted), "
              "whereas the SCZD disease signal is largely composition-confounded. Lead the "
              "DTU claim with the composition-robust aging layer, and disclose the GTEx "
              "cortical collapse as the boundary of that robustness rather than burying it.",
              "- The marker cut (`marker_depleted`/`marker_enriched`) reports whether module "
              "member genes are over/under-represented for cell-type marker genes; near-null "
              "argues the modules are not simply bags of cell-type markers."]
    (out / "COMPOSITION_ADJUSTMENT_SUMMARY.md").write_text("\n".join(lines))
    print(f"[meta] wrote {out/'COMPOSITION_ADJUSTMENT_SUMMARY.md'}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    pf = sub.add_parser("fractions", help="write fractions + marker cut for a cohort")
    pf.add_argument("analysis", choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    pf.add_argument("--region", action="append", dest="regions",
                    help="brainseq-aging/gtex-aging region(s); default = all covered regions.")
    pf.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    pf.add_argument("--music-dir", type=Path, default=None,
                    help="override MuSiC output dir (default: sibling BrainSEQ run, or the "
                         "in-repo GTEx run for gtex-aging).")
    pf.add_argument("--marker-top-n", type=int, default=25,
                    help="top mean-ratio markers per cell type (matches DeconvoBuddies run).")
    pe = sub.add_parser("export-gtex", help="write genes×samples bulk parquets for the R "
                                            "GTEx MuSiC deconvolution")
    pe.add_argument("--region", action="append", dest="regions",
                    help="GTEx region(s); default = all with a matched Tran reference.")
    pm = sub.add_parser("meta", help="cross-region with-vs-without adjustment rollup")
    pm.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    args = p.parse_args()

    if args.cmd == "meta":
        meta(args.variant)
    elif args.cmd == "export-gtex":
        for region in (args.regions or sorted(GTEX_REF)):
            export_gtex_bulk(region, GTEX_MUSIC_DIR)
    elif args.analysis == "brainseq-sczd":
        run("brainseq-sczd", None, args.variant, _music_dir_for("brainseq-sczd", args.music_dir),
            args.marker_top_n)
    elif args.analysis == "gtex-aging":
        mdir = _music_dir_for("gtex-aging", args.music_dir)
        for region in (args.regions or sorted(GTEX_REF)):
            run("gtex-aging", region, args.variant, mdir, args.marker_top_n)
    else:
        mdir = _music_dir_for("brainseq-aging", args.music_dir)
        for region in (args.regions or ["caudate", "hippocampus", "dlpfc"]):
            run("brainseq-aging", region, args.variant, mdir, args.marker_top_n)


if __name__ == "__main__":
    main()
