"""Module-level genetic anchoring — is the co-switch *program* anchored to splicing genetics?

The QTL-anchoring result (`qtl_anchoring.py`) is a *pooled* gene-set test: co-switch genes,
in aggregate, show a splicing-specific sQTL/eQTL contrast. But the manuscript claim is about
genetically anchored *programs* (modules), and a reviewer can fairly object that pooling
across all co-switch genes does not show any individual module is anchored — the coordination
itself is untested (Finding 3 caveat).

This module closes that gap at module granularity, reusing the same power-matched machinery:
for each phenotype-associated module we compute the splicing-specificity contrast on that
module's own member genes (sQTL matched-OR vs eQTL matched-OR, i.e. log-OR difference) against
the method's tested-gene universe, and calibrate it with a **permutation null** of random
equal-size gene sets drawn from the universe. A module whose contrast exceeds the size-matched
null is genetically anchored to splicing *as a unit* — not just carried by the pooled average.

Reuses `qtl_anchoring.{load_qtl_genes, isoform_counts, matched_enrichment, resolve_tissue}`
and the GTEx sGenes/eGenes catalogs already on disk; no individual-level genotypes, so no
controlled-access dependency. Deterministic (seed). This is the tractable module-level test;
the eigenswitch×genotype and module-restricted S-LDSC routes (which need new controlled
genotype extraction) are noted in the summary as deeper follow-ons.

Per analysis writes, under <artifact-parent>/_m/:
  module_genetic_anchoring.parquet — one row per pheno-sig module: sQTL/eQTL OR, contrast,
      permutation p, go_invisible.
  MODULE_GENETIC_ANCHORING.md — the per-module writeup.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.interpret_modules import DEFAULT_GTF_CACHE
from isograph_benchmark.real_data.qtl_anchoring import (
    DEFAULT_XQTL_DIR,
    _bare,
    isoform_counts,
    load_qtl_genes,
    matched_enrichment,
    resolve_tissue,
)
from isograph_benchmark.real_data.partition_provenance import load_enrichment
from isograph_benchmark.real_data.sweep_leiden import _artifact_dir

_ANALYSES = {
    "brainseq-sczd": [None],
    "brainseq-aging": ["caudate", "hippocampus", "dlpfc"],
    "gtex-aging": None,  # resolved to GTEX_REGIONS at call time
}


def _qtl_universe(xqtl_dir: Path, tissue: str, kind: str, fdr: float,
                  iso_genes: set[str], iso_per_gene: pd.Series) -> pd.DataFrame:
    """Tested-gene universe for one xQTL kind with the QTL-detectability covariates
    (mirrors qtl_anchoring.run_anchoring's inner loop)."""
    qtl = load_qtl_genes(xqtl_dir, tissue, kind, fdr)
    qtl = qtl[qtl["gene"].isin(iso_genes)].copy()
    qtl["log_num_var"] = np.log1p(qtl["num_var"])
    qtl["log_len"] = np.log(qtl["gene_len"])
    qtl["log_iso"] = np.log(qtl["gene"].map(iso_per_gene).fillna(1.0).clip(lower=1))
    if "group_size" in qtl.columns:
        qtl["log_group"] = np.log(qtl["group_size"].clip(lower=1))
    return qtl


def _contrast(sqtl: pd.DataFrame, eqtl: pd.DataFrame, genes: set[str]) -> tuple[float, dict, dict]:
    """log(sQTL matched-OR) - log(eQTL matched-OR) for a gene set; NaN if unestimable."""
    s = matched_enrichment(sqtl, genes)
    e = matched_enrichment(eqtl, genes)
    so, eo = s.get("odds_ratio"), e.get("odds_ratio")
    if so is None or eo is None or not (np.isfinite(so) and np.isfinite(eo)) or so <= 0 or eo <= 0:
        return np.nan, s, e
    return float(np.log(so) - np.log(eo)), s, e


def run_module_anchoring(analysis: str, region: str | None, variant: str, fdr: float,
                         xqtl_dir: Path, gtf_cache: Path, n_perm: int, seed: int) -> pd.DataFrame:
    iso_dir = _artifact_dir(analysis, region, variant)
    mod_dir = iso_dir  # isograph_vae
    tissue = resolve_tissue(analysis, region)

    modules = pd.read_parquet(mod_dir / "modules.parquet")
    modules["gene"] = _bare(modules["gene_id"])
    modules["module_id"] = modules["module_id"].astype(str)
    iso_genes = set(modules["gene"])
    iso_per_gene = isoform_counts(gtf_cache)

    enrich_path = iso_dir.parent / "module_enrichment" / "isograph_modules.parquet"
    enrich = load_enrichment(
        enrich_path, modules, context=f"module_genetic_anchoring {analysis}/{region or ''}"
    )
    if enrich is None:
        print(f"[{analysis}/{region}] no module_enrichment table — skipping", flush=True)
        return pd.DataFrame()
    enrich["module_id"] = enrich["module_id"].astype(str)
    sig = enrich[enrich["pheno_fdr"] <= fdr]
    if sig.empty:
        print(f"[{analysis}/{region}] no phenotype-significant modules — skipping", flush=True)
        return pd.DataFrame()

    sqtl = _qtl_universe(xqtl_dir, tissue, "sQTL", fdr, iso_genes, iso_per_gene)
    eqtl = _qtl_universe(xqtl_dir, tissue, "eQTL", fdr, iso_genes, iso_per_gene)
    universe = sorted(iso_genes)
    rng = np.random.default_rng(seed)

    rows = []
    for _, mrow in sig.iterrows():
        mid = mrow["module_id"]
        genes = set(modules.loc[modules["module_id"] == mid, "gene"])
        obs, s, e = _contrast(sqtl, eqtl, genes)
        go_invisible = bool(mrow.get("n_go_terms", 1) == 0)
        rec = {"analysis": analysis, "region": region or "", "tissue": tissue,
               "module_id": mid, "module_size": len(genes), "go_invisible": go_invisible,
               "sqtl_or": s.get("odds_ratio"), "eqtl_or": e.get("odds_ratio"),
               "n_fg_sqtl": s.get("n_foreground"), "n_fg_eqtl": e.get("n_foreground"),
               "contrast_log": obs}
        if np.isfinite(obs) and n_perm > 0:
            null = np.empty(n_perm)
            k = len(genes)
            for i in range(n_perm):
                samp = set(rng.choice(universe, size=k, replace=False))
                null[i], _, _ = _contrast(sqtl, eqtl, samp)
            null = null[np.isfinite(null)]
            rec["perm_p"] = float((np.sum(null >= obs) + 1) / (len(null) + 1)) if len(null) else np.nan
            rec["null_mean"] = float(np.mean(null)) if len(null) else np.nan
            rec["n_perm_valid"] = int(len(null))
        else:
            rec["perm_p"] = np.nan
            rec["null_mean"] = np.nan
            rec["n_perm_valid"] = 0
        rows.append(rec)
        print(f"[{analysis}/{region}] {mid} (n={len(genes)}, "
              f"{'GO-inv' if go_invisible else 'GO-vis'}): contrast={obs:.3f} "
              f"perm_p={rec['perm_p']}", flush=True)

    res = pd.DataFrame(rows)
    out_dir = ensure_dir(iso_dir.parent)
    res.to_parquet(out_dir / "module_genetic_anchoring.parquet", index=False, compression="zstd")
    _write_report(out_dir, analysis, region, tissue, res, n_perm)
    return res


def _md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    head = "| " + " | ".join(map(str, cols)) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |"
            for r in df.itertuples(index=False)]
    return "\n".join([head, sep, *body])


def _write_report(out_dir: Path, analysis: str, region: str | None, tissue: str,
                  res: pd.DataFrame, n_perm: int) -> None:
    show = res.copy()
    for c in ("sqtl_or", "eqtl_or", "contrast_log", "null_mean"):
        if c in show:
            show[c] = show[c].round(3)
    if "perm_p" in show:
        show["perm_p"] = show["perm_p"].apply(lambda p: f"{p:.3f}" if pd.notna(p) else "NA")
    keep = ["module_id", "module_size", "go_invisible", "sqtl_or", "eqtl_or",
            "contrast_log", "null_mean", "perm_p"]
    lines = [
        f"# Module-level genetic anchoring — {analysis}" + (f"/{region}" if region else "")
        + f" vs GTEx {tissue} xQTL", "",
        "Per phenotype-associated module: the splicing-specificity contrast "
        "(log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, "
        f"calibrated against a size-matched permutation null ({n_perm} draws). "
        "A positive contrast with small `perm_p` = the module is anchored to *splicing* "
        "genetics as a unit, beyond a random equal-size gene set.", "",
        "Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring "
        f"--analysis {analysis}" + (f" --region {region}" if region else "") + "`", "",
        _md_table(show[keep]), "",
        "_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the "
        "module concentrates splicing-anchored genes above chance, not that a single variant "
        "drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests "
        "(new controlled-genotype extraction) are the deeper follow-ons._",
    ]
    (out_dir / "MODULE_GENETIC_ANCHORING.md").write_text("\n".join(lines))
    print(f"[{analysis}/{region}] wrote {out_dir/'MODULE_GENETIC_ANCHORING.md'}", flush=True)


def _artifact_parent(analysis: str, region: str | None, variant: str) -> Path:
    return _artifact_dir(analysis, region, variant).parent


def meta(variant: str) -> None:
    """Pool the 17 per-analysis module_genetic_anchoring.parquet files into one table +
    a Stouffer-style headline: are pheno-sig modules anchored to splicing above the null,
    split by GO-invisible vs GO-visible?"""
    from isograph_benchmark.real_data.run_models import GTEX_REGIONS
    specs = ([("brainseq-sczd", None)]
             + [("brainseq-aging", r) for r in ("caudate", "hippocampus", "dlpfc")]
             + [("gtex-aging", r) for r in GTEX_REGIONS])
    frames = []
    for analysis, region in specs:
        p = _artifact_parent(analysis, region, variant) / "module_genetic_anchoring.parquet"
        if p.exists():
            frames.append(pd.read_parquet(p))
    if not frames:
        print("[meta] no per-analysis parquets found", flush=True)
        return
    allm = pd.concat(frames, ignore_index=True)
    out = ensure_dir(stage_out("anchoring", "module_genetic_anchoring_meta"))
    allm.to_parquet(out / "module_genetic_anchoring_all.parquet", index=False, compression="zstd")

    est = allm[np.isfinite(allm["contrast_log"]) & allm["perm_p"].notna()].copy()

    def _blk(df: pd.DataFrame) -> dict:
        n = len(df)
        return {
            "n_modules": n,
            "n_pos_contrast": int((df["contrast_log"] > 0).sum()),
            "n_anchored_p05": int((df["perm_p"] <= 0.05).sum()),
            "median_contrast": round(float(df["contrast_log"].median()), 3) if n else np.nan,
            "frac_pos": round(float((df["contrast_log"] > 0).mean()), 2) if n else np.nan,
        }

    strata = {
        "all_phenosig": est,
        "go_invisible": est[est["go_invisible"] == True],
        "go_visible": est[est["go_invisible"] == False],
        "neurodegen_gtex": est[est["analysis"] == "gtex-aging"],
        "brainseq_sczd": est[est["analysis"] == "brainseq-sczd"],
    }
    tab = pd.DataFrame({k: _blk(v) for k, v in strata.items()}).T.reset_index().rename(
        columns={"index": "stratum"})
    tab.to_parquet(out / "module_genetic_anchoring_meta.parquet", index=False, compression="zstd")

    top = (est.sort_values("perm_p")
           .head(20)[["analysis", "region", "module_id", "module_size", "go_invisible",
                      "sqtl_or", "eqtl_or", "contrast_log", "perm_p"]].round(3))
    lines = [
        "# Module-level genetic anchoring — cross-analysis meta", "",
        "Per phenotype-associated module, the splicing-specificity contrast "
        "(log sQTL matched-OR − log eQTL matched-OR) vs a size-matched permutation null. "
        "Pooled over 17 analyses (SCZD + 3 BrainSEQ aging + 13 GTEx aging).", "",
        "## By stratum", "", _md_table(tab), "",
        "## Most-anchored modules (smallest perm_p)", "", _md_table(top), "",
        "_A module with positive contrast + small perm_p is anchored to splicing genetics as a "
        "unit — evidence the *program*, not just pooled member genes, is genetically anchored. "
        "SCZD is expected to be weak (disease is eQTL-led); the splicing-anchored modules should "
        "concentrate in GO-invisible / neurodegen-relevant GTEx analyses._",
    ]
    (out / "MODULE_GENETIC_ANCHORING_META.md").write_text("\n".join(lines))
    print(f"[meta] {len(est)} estimable modules; "
          f"anchored p<=0.05: {int((est['perm_p']<=0.05).sum())} "
          f"(GO-inv {int((est[est['go_invisible']==True]['perm_p']<=0.05).sum())})", flush=True)
    print(f"[meta] wrote {out/'MODULE_GENETIC_ANCHORING_META.md'}", flush=True)


def main() -> None:
    from isograph_benchmark.real_data.run_models import GTEX_REGIONS
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--meta", action="store_true",
                   help="pool per-analysis parquets into the cross-analysis meta + summary.")
    p.add_argument("--analysis",
                   choices=["brainseq-sczd", "brainseq-aging", "gtex-aging"])
    p.add_argument("--region", action="append", dest="regions")
    p.add_argument("--variant", choices=["standard", "with-abundance"], default="standard")
    p.add_argument("--fdr", type=float, default=0.10)
    p.add_argument("--xqtl-dir", type=Path, default=DEFAULT_XQTL_DIR)
    p.add_argument("--gtf-cache", type=Path, default=DEFAULT_GTF_CACHE)
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=13)
    args = p.parse_args()

    if args.meta:
        meta(args.variant)
        return
    if not args.analysis:
        p.error("--analysis is required unless --meta")
    if args.analysis == "brainseq-sczd":
        regions = [None]
    elif args.analysis == "gtex-aging":
        regions = args.regions or GTEX_REGIONS
    else:
        regions = args.regions or ["caudate", "hippocampus", "dlpfc"]
    for region in regions:
        run_module_anchoring(args.analysis, region, args.variant, args.fdr,
                             args.xqtl_dir, args.gtf_cache, args.n_perm, args.seed)


if __name__ == "__main__":
    main()
