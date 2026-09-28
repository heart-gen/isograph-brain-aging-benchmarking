"""How age-coupled is estimated cell-type composition, and does that explain the loss?

The composition-adjustment contrast (``celltype_composition meta``) shows that adjusting
for MuSiC cell-type fractions removes most switch-unique aging genes in the two GTEx
cortical regions while leaving limbic/striatal regions largely intact. That result is
reported both ways in the manuscript -- composition as confounder, or composition as part
of the aging process -- but nothing in the repository yet *shows* how strongly composition
itself tracks age in each region, so the reader cannot judge which regions are at risk of
over-adjustment.

This module supplies that missing panel. It is a pure read over artifacts already on disk:
the committed per-sample MuSiC fractions (``<artifact_dir>/celltype_fractions.parquet``,
written by ``celltype_composition fractions``) joined to the donor age already carried in
each bundle's sample table. No deconvolution and no model fitting is redone here.

Three outputs, all keyed by the same analysis labels the composition rollup uses:

1. ``composition_age_coupling.parquet`` -- one row per analysis x cell type: Spearman rho
   of estimated proportion against age, its p value and a within-analysis BH FDR.
2. ``composition_age_samples.csv.gz`` -- the tidy per-sample fractions + age behind (1),
   so the figure can draw the scatter without re-deriving the join.
3. ``composition_age_persistence.csv`` -- each analysis's strongest age-coupled cell type
   against the gene-level persistence of its switch-unique set (``n_overlap`` /
   ``comp_unique_base`` from the composition rollup). This is the panel that makes the
   confounder-vs-mediator caveat concrete: the regions that lose the most are the regions
   whose composition moves most with age.
4. ``composition_marker_tests.csv`` -- per analysis, how many module x cell-type marker
   tests were run and how many reached FDR < 0.05 for enrichment and for depletion. This
   is the control that does not depend on the covariate model at all: if the modules were
   collections of cell-type markers, adjustment would be removing the modules' own
   definition rather than a confounder.

Correlation, not causation: a strong rho does not establish that composition drives the
switch signal, and a weak rho does not license a cell-intrinsic reading. The panel bounds
the over-adjustment risk; it does not resolve it.

Run (local, no scheduler): ``bash 03_module_characterization/_h/04b.composition_age_coupling.sh``
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, rel, stage_out

# Per-sample age lives in the bundle sample table, under a cohort-specific column name --
# the same columns run_models.py fits the age spline on.
AGE_COL = {"brainseq": "Age", "gtex": "AGE"}
DONOR_COL = {"brainseq": "BrNum", "gtex": "SUBJID"}

# Anatomical class, not cohort, is what the composition story splits on; reused verbatim
# from the composition figure so the two supplements colour the same regions alike.
CLASS_OF = {
    "SCZD (caudate)": "Disease (SCZD)",
    "aging caudate": "Limbic / striatal",
    "aging hippocampus": "Limbic / striatal",
    "aging DLPFC": "Cortical",
    "GTEx amygdala": "Limbic / striatal",
    "GTEx hippocampus": "Limbic / striatal",
    "GTEx caudate_basal_ganglia": "Limbic / striatal",
    "GTEx putamen_basal_ganglia": "Limbic / striatal",
    "GTEx nucleus_accumbens_basal_ganglia": "Limbic / striatal",
    "GTEx anterior_cingulate_cortex_ba24": "Cortical",
    "GTEx frontal_cortex_ba9": "Cortical",
    "GTEx cortex": "Cortical",
}


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg, inline so this module needs no statsmodels."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _bundle_path(cohort: str, region: str) -> Path:
    """Bundle for a module-discovery store, mirroring sweep_leiden._bundle_path.

    Resolved here rather than imported so this module stays free of the modelling stack.
    """
    if cohort == "brainseq":
        if region == "caudate_sczd":
            return rel("inputs", "bundles", "brainseq_sczd", "caudate")
        return rel("inputs", "bundles", "brainseq_v1", region)
    if cohort == "gtex":
        return rel("inputs", "bundles", "gtex_v11_brain", region)
    raise ValueError(f"Unknown cohort: {cohort!r}")


def _label(cohort: str, region: str) -> str:
    """The analysis label the composition rollup uses, so the two tables join."""
    if cohort == "gtex":
        return f"GTEx {region}"
    if region == "caudate_sczd":
        return "SCZD (caudate)"
    return f"aging {'DLPFC' if region == 'dlpfc' else region}"


def _stores(variant: str) -> list[tuple[str, str, Path]]:
    """Every analysis with committed MuSiC fractions, discovered rather than hard-coded.

    The deconvolved set is defined by what is on disk: regions without a defensibly matched
    snRNA reference were never deconvolved, so they are absent here by construction.
    """
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    found = []
    for p in sorted(stage_out("modules").glob(f"*/*/_m/{subdir}/celltype_fractions.parquet")):
        cohort, region = p.parts[-5], p.parts[-4]
        found.append((cohort, region, p))
    return found


def _samples(cohort: str, region: str, fractions_p: Path) -> pd.DataFrame:
    """Tidy sample x cell-type fractions joined to donor age."""
    fr = pd.read_parquet(fractions_p)
    samples = pd.read_parquet(_bundle_path(cohort, region) / "samples.parquet")
    age_col, donor_col = AGE_COL[cohort], DONOR_COL[cohort]
    keep = ["sample_id", age_col] + ([donor_col] if donor_col in samples.columns else [])
    meta = samples[keep].rename(columns={age_col: "age", donor_col: "donor_id"})
    cell_cols = [c for c in fr.columns if c != "sample_id"]
    long = fr.melt(id_vars="sample_id", value_vars=cell_cols,
                   var_name="cell_type", value_name="proportion")
    out = long.merge(meta, on="sample_id", how="inner")
    out.insert(0, "label", _label(cohort, region))
    out.insert(1, "cohort", "BrainSEQ" if cohort == "brainseq" else "GTEx")
    return out.dropna(subset=["age", "proportion"])


def _coupling(tidy: pd.DataFrame) -> pd.DataFrame:
    """Spearman rho of each cell type's estimated proportion against age.

    Spearman rather than Pearson: MuSiC proportions are bounded, often zero-inflated, and
    the age effect need not be linear on the proportion scale. A cell type that is zero in
    every sample carries no information and is reported with a missing rho rather than
    dropped, so the panel shows the full reference panel that was adjusted for.
    """
    rows = []
    for (label, cell), g in tidy.groupby(["label", "cell_type"], sort=False):
        prop = g["proportion"].to_numpy()
        if len(g) < 3 or np.nanstd(prop) == 0:
            rho, p = np.nan, np.nan
        else:
            rho, p = stats.spearmanr(prop, g["age"].to_numpy())
        rows.append({"label": label, "cohort": g["cohort"].iloc[0], "cell_type": cell,
                     "n_samples": int(len(g)), "mean_proportion": float(np.nanmean(prop)),
                     "rho_age": float(rho), "p_age": float(p)})
    df = pd.DataFrame(rows)
    df["class"] = df["label"].map(CLASS_OF)
    # BH within analysis: each region's reference panel is its own family of tests.
    df["fdr_age"] = np.nan
    for label, g in df.groupby("label"):
        ok = g["p_age"].notna()
        df.loc[g.index[ok], "fdr_age"] = _bh(g.loc[ok, "p_age"].to_numpy())
    return df.sort_values(["label", "cell_type"], ignore_index=True)


def _persistence(coupling: pd.DataFrame) -> pd.DataFrame:
    """Join each analysis's strongest age-coupled cell type to its switch-unique persistence.

    ``n_overlap / comp_unique_base`` is a genuine fraction (the kept genes are a subset of
    the unadjusted set); the ratio of the two counts is not, and is deliberately not used.
    """
    bs = pd.read_parquet(stage_out("characterize") / "composition_adjustment.parquet")
    gt_p = (stage_out("modules", "gtex", "_m", "composition")
            / "composition_adjustment_gtex.parquet")
    parts = [bs] + ([pd.read_parquet(gt_p)] if gt_p.exists() else [])
    roll = pd.concat(parts, ignore_index=True).rename(columns={"region": "label"})

    top = (coupling.dropna(subset=["rho_age"])
           .assign(abs_rho=lambda d: d["rho_age"].abs())
           .sort_values("abs_rho", ascending=False)
           .groupby("label", as_index=False)
           .first()[["label", "cell_type", "rho_age", "abs_rho", "fdr_age"]]
           .rename(columns={"cell_type": "top_cell_type", "rho_age": "top_rho_age",
                            "abs_rho": "max_abs_rho_age", "fdr_age": "top_fdr_age"}))
    # Reference panels differ in size between regions (6-8 cell types), so the count of
    # significant cell types compares granularity as much as age-coupling; carry the
    # fraction alongside it and plot that.
    per_analysis = coupling.groupby("label").agg(
        n_celltypes_tested=("fdr_age", lambda s: int(s.notna().sum())),
        n_celltypes_age_fdr05=("fdr_age", lambda s: int((s < 0.05).sum())),
        median_abs_rho_age=("rho_age", lambda s: float(s.abs().median())),
    ).reset_index()
    per_analysis["frac_celltypes_age_fdr05"] = np.where(
        per_analysis["n_celltypes_tested"] > 0,
        per_analysis["n_celltypes_age_fdr05"] / per_analysis["n_celltypes_tested"], np.nan)

    out = (roll[["label", "comp_unique_base", "comp_unique_adj", "n_overlap", "n_new"]]
           .merge(top, on="label", how="inner").merge(per_analysis, on="label", how="left"))
    out["persistence"] = np.where(out["comp_unique_base"] > 0,
                                  out["n_overlap"] / out["comp_unique_base"], np.nan)
    out["class"] = out["label"].map(CLASS_OF)
    return out.sort_values("max_abs_rho_age", ascending=False, ignore_index=True)


def _marker_tests(variant: str, root: Path | None = None) -> pd.DataFrame:
    """Per-analysis module x cell-type marker tests: how many, and how many significant.

    Read from each store's ``celltype_composition/marker_enrichment.parquet`` (one row per
    module x cell type, with BH-adjusted enrichment and depletion p values). The number of
    tests is carried with the two counts because it differs ~3x between regions -- with the
    module count and the size of the reference panel -- so a bare count of significant
    tests would compare granularity rather than marker content.
    """
    subdir = "isograph_vae_with_abundance" if variant == "with-abundance" else "isograph_vae"
    root = stage_out("modules") if root is None else root
    rows = []
    for f in sorted(root.glob(f"*/*/_m/{subdir}/celltype_composition/"
                              "marker_enrichment.parquet")):
        cohort, region = f.parts[-6], f.parts[-5]
        mk = pd.read_parquet(f, columns=["fdr_enrich", "fdr_deplete"])
        label = _label(cohort, region)
        rows.append(dict(
            label=label, cohort="GTEx" if cohort == "gtex" else "BrainSEQ",
            **{"class": CLASS_OF.get(label)},
            n_marker_tests=int(len(mk)),
            marker_enriched=int((mk["fdr_enrich"] < 0.05).sum()),
            marker_depleted=int((mk["fdr_deplete"] < 0.05).sum()),
        ))
    cols = ["label", "cohort", "class", "n_marker_tests", "marker_enriched",
            "marker_depleted"]
    return pd.DataFrame(rows, columns=cols)


def _md_table(df: pd.DataFrame) -> str:
    head = "| " + " | ".join(df.columns) + " |"
    rule = "| " + " | ".join("---" for _ in df.columns) + " |"
    body = ["| " + " | ".join("" if pd.isna(v) else
                              (f"{v:.3g}" if isinstance(v, float) else str(v))
                              for v in row) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, rule, *body])


def run(variant: str) -> None:
    stores = _stores(variant)
    if not stores:
        raise SystemExit("no celltype_fractions.parquet found; run `celltype_composition "
                         "fractions` first")
    tidy = pd.concat([_samples(c, r, p) for c, r, p in stores], ignore_index=True)
    coupling = _coupling(tidy)
    persist = _persistence(coupling)

    out = ensure_dir(stage_out("characterize"))
    coupling.to_parquet(out / "composition_age_coupling.parquet", index=False,
                        compression="zstd")
    # CSV alongside the parquet: the figure script reads these, so the figure can be built
    # with an `arrow` that was compiled without zstd support.
    coupling.to_csv(out / "composition_age_coupling.csv", index=False)
    # gzipped: this is the one large table (one row per sample x cell type) and R reads
    # .csv.gz directly, so there is no reason to keep it uncompressed in the repo.
    tidy.to_csv(out / "composition_age_samples.csv.gz", index=False, compression="gzip")
    persist.to_csv(out / "composition_age_persistence.csv", index=False)
    markers = _marker_tests(variant)
    markers.to_csv(out / "composition_marker_tests.csv", index=False)

    n_sig = int((coupling["fdr_age"] < 0.05).sum())
    cortical = persist[persist["class"] == "Cortical"]["max_abs_rho_age"]
    limbic = persist[persist["class"] == "Limbic / striatal"]["max_abs_rho_age"]
    summary = {
        "variant": variant, "n_analyses": int(coupling["label"].nunique()),
        "n_celltype_tests": int(len(coupling)),
        "n_celltype_age_fdr05": n_sig,
        "median_max_abs_rho_cortical": float(cortical.median()) if len(cortical) else None,
        "median_max_abs_rho_limbic": float(limbic.median()) if len(limbic) else None,
        "n_marker_tests": int(markers["n_marker_tests"].sum()),
        "n_marker_enriched": int(markers["marker_enriched"].sum()),
        "n_marker_depleted": int(markers["marker_depleted"].sum()),
    }
    (out / "composition_age_coupling_summary.json").write_text(json.dumps(summary, indent=2))

    lines = [
        "# Is estimated cell-type composition itself age-coupled?", "",
        "Spearman correlation of each MuSiC cell-type proportion with donor age, per "
        "analysis, read from the committed fractions and the bundle sample tables. Nothing "
        "is re-deconvolved and no model is refit.", "",
        f"{n_sig} of {len(coupling)} analysis x cell-type tests reach FDR < 0.05 "
        f"(Benjamini-Hochberg within analysis) across {coupling['label'].nunique()} "
        "deconvolved analyses.", "",
        "## Strongest age-coupled cell type vs switch-unique persistence", "",
        _md_table(persist[["label", "class", "top_cell_type", "top_rho_age",
                           "median_abs_rho_age", "n_celltypes_age_fdr05",
                           "n_celltypes_tested", "comp_unique_base", "n_overlap",
                           "persistence"]]), "",
        "## Interpretation", "",
        "- `persistence` is `n_overlap / comp_unique_base`: the share of a region's "
        "unadjusted switch-unique genes that are still switch-unique after adjustment. It "
        "is a real fraction; the ratio of the two counts is not, because adjustment adds "
        "genes as well as dropping them.",
        "- **A strong `top_rho_age` marks a region where adjustment and the age term "
        "compete for the same variance.** In those regions a collapse in switch-unique "
        "genes is equally consistent with composition confounding the signal and with "
        "over-adjustment removing real signal; this table bounds the risk rather than "
        "resolving it.",
        "- **Age-coupling alone does not predict the collapse.** Cortical regions do carry "
        "the most age-coupled reference panels, but ACC BA24 and BrainSEQ DLPFC are "
        "cortical and retain a third of their switch-unique genes, while GTEx BA9 and "
        "cortex retain almost none. Composition-age coupling is a necessary part of the "
        "over-adjustment story, not a sufficient one; the reference match matters too "
        "(GTEx cortex/BA9 are deconvolved against a DLPFC snRNA panel).",
        "- Correlation is not causation in either direction. Read this panel next to the "
        "marker cut below, which does not depend on the covariate model at all.",
        "",
        "## Are the modules collections of cell-type markers?", "",
        f"{int(markers['marker_enriched'].sum())} of {int(markers['n_marker_tests'].sum())} "
        "module x cell-type tests reach FDR < 0.05 for marker enrichment and "
        f"{int(markers['marker_depleted'].sum())} for depletion "
        "(`composition_marker_tests.csv`).", "",
        _md_table(markers),
    ]
    (out / "COMPOSITION_AGE_COUPLING.md").write_text("\n".join(lines))
    print(f"[coupling] {coupling['label'].nunique()} analyses, {n_sig}/{len(coupling)} "
          f"cell-type tests at FDR<0.05", flush=True)
    print(f"[coupling] wrote {out/'COMPOSITION_AGE_COUPLING.md'}", flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    run(p.parse_args().variant)


if __name__ == "__main__":
    main()
