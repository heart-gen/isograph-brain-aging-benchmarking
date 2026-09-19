"""Do IsoGraph's switch genes and switched exons carry clinical/selective-constraint signal?

The switch coding-consequence analysis showed the switch axis is productive UTR/CDS
remodeling. This adds the orthogonal functional-impact question a reviewer will ask, with two
light-scope, download-only evidence sources:

  1. gnomAD LOEUF (PRIMARY, gene-level anchor). Are the switch genes under stronger loss-of-
     function constraint than genes genome-wide? Mann-Whitney of switch-gene LOEUF vs all genes
     in the constraint table ('less' = switch genes more constrained), stratified GO-invisible
     vs GO-visible. This is the clean gene-level anchor.

  2. ClinVar pathogenic density (SECONDARY, exon-level, within-gene, DIRECTION-NEUTRAL). For each
     switch gene we split its switch-pair exons into SWITCHED (present in exactly one isoform of
     some pair — differentially used) and BACKGROUND (constitutive across the switching isoforms),
     and count Pathogenic / Likely_pathogenic ClinVar variants per kb in each. The test is a
     WITHIN-GENE label permutation (shuffle switched/background labels among a gene's own exons,
     keeping the switched count) so gene-level ClinVar ascertainment cancels exactly — same logic
     as switch_consequence's within-gene null. Effect = switched/background P/LP density ratio;
     p = TWO-sided empirical permutation p. NB alternatively-spliced exons are typically LESS
     constrained than constitutive coding exons, so ratio < 1 (depletion) is the expected null
     biology here, not a defect — the test asks whether IsoGraph's switches deviate from it.

One region+resolution per invocation (mirrors switch_consequence). Deterministic given --seed.
Requires the downloaded ClinVar VCF + gnomAD constraint table (see 06_switch_mechanism/_h/01f.download_clinical.sh).

Output (<artifact_dir>/clinical_consequence/): exon_clinvar.parquet (per exon), gene_constraint
.parquet (per gene LOEUF), clinical_consequence.parquet (per stratum: densities, ratio, perm p),
constraint_summary.parquet (per stratum: median LOEUF, MWU p), CLINICAL_CONSEQUENCE.md.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from isograph.explain.structure import parse_gtf
from isograph_benchmark.paths import ensure_dir, region_store, rel
from isograph_benchmark.real_data.coloc_prep import load_switch_genes
from isograph_benchmark.real_data.interpret_modules import (
    DEFAULT_GTF_CACHE,
    DEFAULT_GTF_PATH,
)
from isograph_benchmark.real_data.qtl_anchoring import _bare

_CLINICAL_DIR = rel("inputs", "raw", "clinical")
_CLINVAR_VCF = _CLINICAL_DIR / "clinvar.vcf.gz"
_GNOMAD_CONSTRAINT = _CLINICAL_DIR / "gnomad.v4.1.constraint_metrics.tsv"
# ClinVar CLNSIG values counted as pathogenic.
_PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
        "Pathogenic,_low_penetrance", "Likely_pathogenic,_low_penetrance"}


def _strip_version(x: str) -> str:
    return x.split(".", 1)[0]


def _norm_chrom(chrom: str) -> str:
    """GTF 'chr1' -> ClinVar '1' (and chrM -> MT)."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return "MT" if c in ("M", "MT") else c


def _load_clinvar() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Parse the ClinVar VCF once into per-chromosome sorted POS arrays: all variants and the
    Pathogenic/Likely_pathogenic subset (for the P/LP density and its all-ClinVar denominator)."""
    if not _CLINVAR_VCF.exists():
        raise SystemExit(f"{_CLINVAR_VCF} missing; run 06_switch_mechanism/_h/01f.download_clinical.sh first.")
    all_pos: dict[str, list[int]] = {}
    plp_pos: dict[str, list[int]] = {}
    with gzip.open(_CLINVAR_VCF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            chrom, pos, info = f[0], f[1], f[7] if len(f) > 7 else ""
            try:
                p = int(pos)
            except ValueError:
                continue
            all_pos.setdefault(chrom, []).append(p)
            clnsig = ""
            for kv in info.split(";"):
                if kv.startswith("CLNSIG="):
                    clnsig = kv[7:]
                    break
            if clnsig in _PLP:
                plp_pos.setdefault(chrom, []).append(p)
    to_arr = lambda d: {c: np.array(sorted(v), dtype=np.int64) for c, v in d.items()}
    return to_arr(all_pos), to_arr(plp_pos)


def _count_in(arr: np.ndarray | None, start: int, end: int) -> int:
    if arr is None or arr.size == 0:
        return 0
    lo = np.searchsorted(arr, start, side="left")
    hi = np.searchsorted(arr, end, side="right")
    return int(hi - lo)


def _merge_intervals(ivs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sort + merge overlapping/adjacent (start, end) intervals."""
    out: list[tuple[int, int]] = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _overlaps(exon: tuple[int, int], merged: list[tuple[int, int]]) -> bool:
    """Does an exon interval intersect a merged (sorted) CDS footprint?"""
    s, e = exon
    for cs, ce in merged:
        if cs > e:
            break
        if s <= ce:
            return True
    return False


def _switch_exon_sets(tx_db: dict, pairs: pd.DataFrame):
    """For one gene's switch pairs: SWITCHED exons (in exactly one isoform of some pair) and
    BACKGROUND exons (all other exons of the switching isoforms). Returns (switched, background,
    chrom, cds_footprint) — exon sets of (start, end) and the merged CDS interval list across the
    switching isoforms (for the coding-exon restriction)."""
    ex_by_tx: dict[str, set] = {}
    cds: list[tuple[int, int]] = []
    chrom = None
    for tx in set(pairs["transcript_id_1"]) | set(pairs["transcript_id_2"]):
        rec = tx_db.get(tx)
        if rec is None:
            continue
        ex_by_tx[tx] = {tuple(e) for e in rec.exons}
        cds.extend(tuple(c) for c in rec.cds)
        chrom = rec.chrom
    switched: set = set()
    for r in pairs.itertuples():
        e1, e2 = ex_by_tx.get(r.transcript_id_1), ex_by_tx.get(r.transcript_id_2)
        if e1 and e2:
            switched |= (e1 ^ e2)          # differentially used exons
    allex: set = set().union(*ex_by_tx.values()) if ex_by_tx else set()
    background = allex - switched
    return switched, background, chrom, _merge_intervals(cds)


def _exon_table(tx_db: dict, obs: pd.DataFrame, plp: dict, allcv: dict) -> pd.DataFrame:
    """Per-exon rows (gene, chrom, start, end, length, switched, cds_overlap, n_plp, n_clinvar,
    go_invisible). `cds_overlap` = exon intersects the gene's CDS footprint (coding exon)."""
    rows = []
    for gene, sub in obs.groupby("gene"):
        go_inv = bool(sub["go_invisible"].iloc[0])
        switched, background, chrom, cds = _switch_exon_sets(tx_db, sub)
        if chrom is None:
            continue
        c = _norm_chrom(chrom)
        for is_sw, exons in ((True, switched), (False, background)):
            for (s, e) in exons:
                rows.append((gene, chrom, int(s), int(e), int(e - s + 1), is_sw,
                             _overlaps((s, e), cds),
                             _count_in(plp.get(c), s, e), _count_in(allcv.get(c), s, e),
                             go_inv))
    return pd.DataFrame(rows, columns=["gene", "chrom", "start", "end", "length", "switched",
                                       "cds_overlap", "n_plp", "n_clinvar", "go_invisible"])


def _density(sub: pd.DataFrame, mask: np.ndarray) -> float:
    kb = sub.loc[mask, "length"].sum() / 1000.0
    return float(sub.loc[mask, "n_plp"].sum() / kb) if kb > 0 else np.nan


def _within_gene_perm(exons: pd.DataFrame, n_perm: int, rng: np.random.Generator) -> dict:
    """Pooled switched vs background P/LP density + one-sided within-gene label-permutation p.
    Only genes with >=1 switched AND >=1 background exon contribute permutable label freedom."""
    sw = exons["switched"].to_numpy()
    obs_ratio_num = _density(exons, sw)
    obs_ratio_den = _density(exons, ~sw)
    ratio = (obs_ratio_num / obs_ratio_den
             if obs_ratio_den and obs_ratio_den > 0 else np.nan)
    # genes eligible for label shuffling
    perm_genes = [g for g, s in exons.groupby("gene")
                  if s["switched"].any() and (~s["switched"]).any()]
    if not perm_genes or np.isnan(ratio):
        return {"switched_plp_per_kb": obs_ratio_num, "bg_plp_per_kb": obs_ratio_den,
                "ratio": ratio, "p_emp": np.nan, "n_perm_genes": len(perm_genes)}
    pe = exons[exons["gene"].isin(perm_genes)].copy()
    length = pe["length"].to_numpy(float)
    nplp = pe["n_plp"].to_numpy(float)
    gene_idx = {g: pe.index[pe["gene"] == g].to_numpy() for g in perm_genes}
    nsw = {g: int(pe.loc[idx, "switched"].sum()) for g, idx in gene_idx.items()}
    obs_sw_density = _density(pe, pe["switched"].to_numpy())
    ge = le = 0
    pos = {v: i for i, v in enumerate(pe.index)}
    for _ in range(n_perm):
        lab = np.zeros(len(pe), dtype=bool)
        for g, idx in gene_idx.items():
            chosen = rng.choice(len(idx), size=nsw[g], replace=False)
            for c in chosen:
                lab[pos[idx[c]]] = True
        kb = length[lab].sum() / 1000.0
        d = nplp[lab].sum() / kb if kb > 0 else 0.0
        ge += int(d >= obs_sw_density)
        le += int(d <= obs_sw_density)
    # two-sided: switched exons may be enriched (ratio>1) OR depleted (ratio<1) for pathogenic
    # variants relative to constitutive exons; alt-spliced exons are typically less constrained.
    p_two = min(1.0, 2.0 * min(ge + 1, le + 1) / (n_perm + 1))
    return {"switched_plp_per_kb": obs_ratio_num, "bg_plp_per_kb": obs_ratio_den,
            "ratio": ratio, "p_emp": p_two, "n_perm_genes": len(perm_genes)}


def _load_constraint() -> pd.DataFrame:
    """gnomAD constraint table -> bare gene_id + LOEUF + missense o/e (one row per gene, the
    most-constrained transcript)."""
    if not _GNOMAD_CONSTRAINT.exists():
        raise SystemExit(f"{_GNOMAD_CONSTRAINT} missing; run 06_switch_mechanism/_h/01f.download_clinical.sh first.")
    df = pd.read_csv(_GNOMAD_CONSTRAINT, sep="\t", dtype=str, low_memory=False)
    gid = next((c for c in ("gene_id", "gene_id.", "ensembl_gene_id", "gene")
                if c in df.columns), None)
    loeuf = next((c for c in ("lof.oe_ci.upper", "oe_lof_upper", "loeuf")
                  if c in df.columns), None)
    mis = next((c for c in ("mis.oe", "oe_mis") if c in df.columns), None)
    if gid is None or loeuf is None:
        raise SystemExit(f"constraint table missing gene/LOEUF columns; has {list(df.columns)[:8]}")
    out = pd.DataFrame({"gene": df[gid].map(_strip_version),
                        "loeuf": pd.to_numeric(df[loeuf], errors="coerce"),
                        "mis_oe": pd.to_numeric(df[mis], errors="coerce") if mis else np.nan})
    out = out.dropna(subset=["loeuf"]).sort_values("loeuf")
    return out.groupby("gene", as_index=False).first()   # most-constrained transcript per gene


def _observed_pairs(artifact_dir: Path, switch_genes: pd.DataFrame) -> pd.DataFrame:
    sp = pd.read_parquet(artifact_dir / "module_interpret" / "structure_switch_pairs.parquet")
    sp["gene"] = _bare(sp["gene_id"])
    tag = switch_genes.groupby("gene")["go_invisible"].max().reset_index()
    sp = sp.merge(tag, on="gene", how="inner")
    return sp[["gene", "transcript_id_1", "transcript_id_2", "go_invisible"]]


def _strata(df: pd.DataFrame, gene_col: str = "gene"):
    yield "all", df
    yield "go_invisible", df[df["go_invisible"]]
    yield "go_visible", df[~df["go_invisible"]]


def run(region: str, artifact_dir: Path, fdr: float, n_perm: int, seed: int) -> pd.DataFrame:
    switch_genes = load_switch_genes(artifact_dir, fdr)
    if switch_genes.empty:
        print(f"{region}: no phenotype-significant switch genes; skipping.")
        return pd.DataFrame()
    obs = _observed_pairs(artifact_dir, switch_genes)
    if obs.empty:
        print(f"{region}: no switch pairs for phenotype-significant genes; skipping.")
        return pd.DataFrame()

    tx_db = parse_gtf(DEFAULT_GTF_PATH, cache=DEFAULT_GTF_CACHE)
    allcv, plp = _load_clinvar()
    exons = _exon_table(tx_db, obs, plp, allcv)
    if exons.empty:
        print(f"{region}: no exons resolvable; skipping.")
        return pd.DataFrame()

    rng = np.random.default_rng(seed)
    rows = []
    # `all_exons` = every switch-pair exon; `cds` = coding exons only (both switched and
    # background restricted to exons overlapping the gene's CDS footprint), which removes the
    # UTR/alt-exon length bias so the contrast is coding-vs-coding.
    for scope, escope in (("all_exons", exons), ("cds", exons[exons["cds_overlap"]])):
        for stratum, sub in _strata(escope):
            if sub.empty:
                continue
            r = _within_gene_perm(sub, n_perm, rng)
            r.update(region=region, scope=scope, stratum=stratum,
                     n_genes=int(sub["gene"].nunique()),
                     n_switched_exons=int(sub["switched"].sum()),
                     n_bg_exons=int((~sub["switched"]).sum()),
                     n_plp_switched=int(sub.loc[sub["switched"], "n_plp"].sum()),
                     n_plp_bg=int(sub.loc[~sub["switched"], "n_plp"].sum()))
            rows.append(r)
    summary = pd.DataFrame(rows)

    # gnomAD LOEUF anchor
    constraint = _load_constraint()
    gene_tag = obs.drop_duplicates("gene")[["gene", "go_invisible"]]
    gc = gene_tag.merge(constraint, on="gene", how="left")
    gc["region"] = region
    all_loeuf = constraint["loeuf"].to_numpy()
    crows = []
    for stratum, sub in _strata(gc):
        vals = sub["loeuf"].dropna().to_numpy()
        mwu_p = (float(mannwhitneyu(vals, all_loeuf, alternative="less").pvalue)
                 if vals.size >= 3 else np.nan)   # 'less' = switch genes MORE constrained
        crows.append({"region": region, "stratum": stratum,
                      "n_genes_loeuf": int(vals.size),
                      "median_loeuf_switch": float(np.median(vals)) if vals.size else np.nan,
                      "median_loeuf_all": float(np.median(all_loeuf)),
                      "mwu_p_more_constrained": mwu_p})
    constraint_summary = pd.DataFrame(crows)

    out_dir = ensure_dir(artifact_dir / "clinical_consequence")
    exons.to_parquet(out_dir / "exon_clinvar.parquet", index=False)
    gc.to_parquet(out_dir / "gene_constraint.parquet", index=False)
    summary.to_parquet(out_dir / "clinical_consequence.parquet", index=False)
    constraint_summary.to_parquet(out_dir / "constraint_summary.parquet", index=False)
    _write_report(out_dir, region, summary, constraint_summary, exons)

    for scope in ("all_exons", "cds"):
        a = summary[(summary["stratum"] == "all") & (summary["scope"] == scope)]
        if not a.empty:
            r = a.iloc[0]
            print(f"{region} [{scope}]: P/LP switched {r['switched_plp_per_kb']:.3f} vs bg "
                  f"{r['bg_plp_per_kb']:.3f} /kb, ratio {r['ratio']:.2f}, p={r['p_emp']:.3g}")
    return summary


def _write_report(out_dir: Path, region: str, summary: pd.DataFrame,
                  constraint: pd.DataFrame, exons: pd.DataFrame) -> None:
    lines = [
        f"# Clinical-consequence of switched exons — {region}", "",
        "**ClinVar P/LP density** (Pathogenic / Likely_pathogenic variants per kb) in exons "
        "that are SWITCHED (differentially used between a switch pair's isoforms) vs BACKGROUND "
        "(constitutive exons of the switching isoforms); `ratio` = switched / background; "
        "`p_emp` = TWO-sided within-gene label-permutation p (gene-level ClinVar ascertainment "
        "cancels). Alt-spliced exons are usually less constrained, so ratio < 1 is the expected "
        "baseline; the gnomAD LOEUF panel below is the primary gene-level anchor. Stratified by "
        "GO-invisible module membership.", "",
        f"- exons scored: **{len(exons)}** "
        f"(switched {int(exons['switched'].sum())}, background {int((~exons['switched']).sum())}; "
        f"CDS-overlapping {int(exons['cds_overlap'].sum())}). `scope` = all switch-pair exons vs "
        f"coding (CDS-overlapping) exons only.",
        "",
        "| scope | stratum | genes | switched/kb | bg/kb | ratio | p | perm genes |",
        "|-------|---------|-------|-------------|-------|-------|---|------------|",
    ]
    for r in summary.itertuples():
        lines.append(f"| {r.scope} | {r.stratum} | {r.n_genes} | {r.switched_plp_per_kb:.3f} | "
                     f"{r.bg_plp_per_kb:.3f} | {r.ratio:.2f} | {r.p_emp:.3g} | {r.n_perm_genes} |")
    lines += ["", "**gnomAD LOEUF** of the switch genes vs all genes in the constraint table "
              "(Mann-Whitney, 'more constrained' = lower LOEUF). Gene-level anchor.", "",
              "| stratum | genes w/ LOEUF | median LOEUF (switch) | median (all) | MWU p |",
              "|---------|----------------|-----------------------|--------------|-------|"]
    for r in constraint.itertuples():
        lines.append(f"| {r.stratum} | {r.n_genes_loeuf} | {r.median_loeuf_switch:.3f} | "
                     f"{r.median_loeuf_all:.3f} | {r.mwu_p_more_constrained:.3g} |")
    (out_dir / "CLINICAL_CONSEQUENCE.md").write_text("\n".join(lines) + "\n")


def _default_artifact_dir(region: str, resolution: str) -> Path:
    sub = "brainseq" if region in {"caudate", "caudate_sczd", "hippocampus", "dlpfc"} else "gtex"
    return region_store(sub, region, resolution)


def main() -> None:
    p = argparse.ArgumentParser(
        description="Within-gene ClinVar/gnomAD clinical-consequence test for switched exons.")
    p.add_argument("--region", required=True, help="tree-qualified region label, e.g. gtex_cortex.")
    p.add_argument("--resolution", default="isograph_vae")
    p.add_argument("--artifact-dir", default=None)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=13)
    args = p.parse_args()
    artifact_dir = (Path(args.artifact_dir) if args.artifact_dir
                    else _default_artifact_dir(args.region, args.resolution))
    if not (artifact_dir / "module_interpret" / "structure_switch_pairs.parquet").exists():
        raise SystemExit(f"no structure_switch_pairs under {artifact_dir}; run interpret_modules first.")
    run(args.region, artifact_dir, args.fdr, args.n_perm, args.seed)


if __name__ == "__main__":
    main()
