"""Three-baseline synthesis: IsoGraph vs the matched WGCNA baselines.

Assembles the head-to-head comparison the matched-feature baselines were built to
settle: does handing WGCNA the *same* switch/multiplex feature matrix IsoGraph
consumes close IsoGraph's gap, or is the network inference (VAE + Leiden) doing
real work beyond the input representation? Four methods per region:
  isograph        — switch+abundance features, VAE + Leiden
  wgcna_gene      — abundance features, classical WGCNA (the input control)
  wgcna_switch_only — switch-only features, WGCNA (matched features, classic infer)
  wgcna_multiplex — switch+abundance features, WGCNA (matched features)

Per region/method it reports module count, median size, GO-enrichment fraction,
phenotype-significant count, and the BOTH count (GO-enriched AND phenotype-sig).
Pools across regions per method, and pulls the cross-cohort GO replication
consistency (isograph vs classical wgcna_gene only — the matched baselines were
not run through replication).

Reads, per region, <_m>/module_enrichment/{<method>_modules.parquet, summary.json}
and 03_module_trust/_m/replication/replication_go_summary.parquet.

Writes under real_data/_m/baseline_comparison/:
  baseline_comparison.parquet        — one row per (cohort, region, method).
  baseline_comparison_pooled.parquet — one row per method, pooled across regions.
  BASELINE_COMPARISON.md             — comparison table + honest win/tie/loss read.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import cohort_dir, ensure_dir, stage_out

# canonical label -> per-region module table basename
METHOD_FILES = {
    "isograph": "isograph_modules.parquet",
    "wgcna_gene": "wgcna_modules.parquet",
    "wgcna_switch_only": "wgcna_switch_modules.parquet",
    "wgcna_multiplex": "wgcna_multiplex_modules.parquet",
}
METHOD_ORDER = list(METHOD_FILES)
# input representation each method is fed (for the features-vs-method read)
FEATURES = {
    "isograph": "switch+abundance",
    "wgcna_gene": "abundance",
    "wgcna_switch_only": "switch-only",
    "wgcna_multiplex": "switch+abundance",
}
# replication_go_summary method labels -> canonical
REPL_METHOD = {"isograph_vae": "isograph", "wgcna_gene": "wgcna_gene"}


def _regions(cohorts: tuple[str, ...]) -> list[tuple[str, str, Path]]:
    out = []
    for cohort in cohorts:
        base = cohort_dir(cohort)
        if not base.is_dir():
            continue
        for region in sorted(d.name for d in base.iterdir() if d.is_dir()):
            enrich = base / region / "_m" / "module_enrichment"
            if (enrich / "summary.json").exists():
                out.append((cohort, region, enrich))
    return out


def _method_row(cohort: str, region: str, method: str, enrich: Path,
                fdr: float) -> dict | None:
    path = enrich / METHOD_FILES[method]
    if not path.exists():
        return None
    m = pd.read_parquet(path)
    pheno = pd.to_numeric(m["pheno_fdr"], errors="coerce")
    go = pd.to_numeric(m["n_go_terms"], errors="coerce").fillna(0)
    sig = pheno <= fdr
    enr = go > 0
    n = len(m)
    n_sig, n_both = int(sig.sum()), int((sig & enr).sum())
    return {
        "cohort": cohort, "region": region, "method": method,
        "features": FEATURES[method], "n_modules": n,
        "median_module_size": float(pd.to_numeric(m["n_genes"]).median()),
        "n_go_enriched": int(enr.sum()),
        "frac_go_enriched": round(float(enr.mean()), 3) if n else np.nan,
        "n_pheno_sig": n_sig, "n_both": n_both,
        "frac_pheno_sig": round(n_sig / n, 3) if n else np.nan,
        "frac_both": round(n_both / n, 3) if n else np.nan,
    }


def collect(cohorts: tuple[str, ...], fdr: float) -> pd.DataFrame:
    rows = []
    for cohort, region, enrich in _regions(cohorts):
        for method in METHOD_ORDER:
            row = _method_row(cohort, region, method, enrich, fdr)
            if row is not None:
                rows.append(row)
    return pd.DataFrame(rows)


def _pool(per: pd.DataFrame) -> pd.DataFrame:
    g = per.groupby("method", sort=False)
    pooled = g.agg(
        n_regions=("region", "size"),
        median_n_modules=("n_modules", "median"),
        median_module_size=("median_module_size", "median"),
        mean_frac_pheno_sig=("frac_pheno_sig", "mean"),
        mean_frac_both=("frac_both", "mean"),
        mean_frac_go_enriched=("frac_go_enriched", "mean"),
        total_pheno_sig=("n_pheno_sig", "sum"),
        total_both=("n_both", "sum"),
    ).reset_index()
    pooled["features"] = pooled["method"].map(FEATURES)
    for c in ("mean_frac_pheno_sig", "mean_frac_both", "mean_frac_go_enriched"):
        pooled[c] = pooled[c].round(3)
    pooled["median_module_size"] = pooled["median_module_size"].round(0)
    order = {m: i for i, m in enumerate(METHOD_ORDER)}
    return pooled.sort_values("method", key=lambda s: s.map(order))


def _replication() -> pd.DataFrame:
    path = stage_out("trust.replication", "replication_go_summary.parquet")
    if not path.exists():
        return pd.DataFrame()
    r = pd.read_parquet(path)
    r["method"] = r["method"].map(REPL_METHOD).fillna(r["method"])
    return r


def run(cohorts: tuple[str, ...] = ("brainseq", "gtex"), fdr: float = 0.10) -> pd.DataFrame:
    per = collect(cohorts, fdr)
    out_dir = ensure_dir(stage_out("characterize", "baseline_comparison"))
    if per.empty:
        print("no module_enrichment summaries found")
        return per
    per.to_parquet(out_dir / "baseline_comparison.parquet", index=False, compression="zstd")
    pooled = _pool(per)
    pooled.to_parquet(out_dir / "baseline_comparison_pooled.parquet", index=False,
                      compression="zstd")
    _write_report(out_dir, per, pooled, _replication(), fdr)
    return pooled


def _markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    rows += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(rows)


def _write_report(out_dir: Path, per: pd.DataFrame, pooled: pd.DataFrame,
                  repl: pd.DataFrame, fdr: float) -> None:
    n_regions = per[["cohort", "region"]].drop_duplicates().shape[0]
    n_full = (per.groupby(["cohort", "region"])["method"].nunique() == len(METHOD_ORDER)).sum()
    ptab = pooled.rename(columns={"mean_frac_pheno_sig": "pheno_sig_rate",
                                  "mean_frac_both": "both_rate",
                                  "mean_frac_go_enriched": "go_rate",
                                  "median_n_modules": "med_n_mod",
                                  "median_module_size": "med_size",
                                  "total_pheno_sig": "tot_pheno_sig"})
    ptab = ptab[["method", "features", "n_regions", "med_n_mod", "med_size",
                 "pheno_sig_rate", "both_rate", "go_rate", "tot_pheno_sig"]]
    cau = per[(per.cohort == "brainseq") & (per.region == "caudate")][
        ["method", "features", "n_modules", "frac_pheno_sig", "frac_both",
         "frac_go_enriched"]]

    lines = [
        "# Three-baseline comparison — IsoGraph vs matched WGCNA baselines",
        "",
        f"Per-region module_enrichment across {n_regions} analyses ({n_full} with all "
        f"four methods). Phenotype-significant = `pheno_fdr <= {fdr}`; GO-enriched = "
        "`n_go_terms > 0`; BOTH = phenotype-significant AND GO-enriched. Rates are "
        "per-module fractions averaged across regions; raw totals scale with module "
        "count (IsoGraph runs at finer resolution) and are NOT directly comparable.",
        "",
        "Reproduce: `python -m isograph_benchmark.real_data.baseline_comparison`.",
        "",
        "## Pooled across regions (per method)",
        "",
        _markdown_table(ptab),
        "",
        "## Caudate (representative single region)",
        "",
        _markdown_table(cau),
        "",
        "## The features-vs-method question",
        "",
        "`wgcna_switch_only` and `wgcna_multiplex` are fed the SAME switch / "
        "switch+abundance feature matrix as IsoGraph but inferred with classical WGCNA, "
        "so the comparison isolates network inference (VAE + Leiden) from the input "
        "representation. `wgcna_gene` is the abundance-only input control.",
        "",
        "## Honest read (per-module rates, not totals)",
        "",
        "- **The phenotype signal lives in the switch features, not the method.** Both "
        "switch-fed methods (isograph, wgcna_switch_only) carry a higher "
        "phenotype-significant rate than the abundance-fed ones (wgcna_gene, "
        "wgcna_multiplex). Representing isoform switching is what buys phenotype "
        "sensitivity.",
        "- **IsoGraph is NOT globally superior on module-level metrics.** Classical "
        "`wgcna_switch_only` matches or exceeds IsoGraph's phenotype-significant rate, "
        "and IsoGraph has the LOWEST BOTH rate (few modules are both phenotype-sig and "
        "GO-enriched, because its GO-enrichment is low by construction). This is the "
        "expected picture: abundance dominates module-level enrichment.",
        "- **One clean method effect survives:** on IDENTICAL switch+abundance features, "
        "IsoGraph's phenotype-significant rate exceeds `wgcna_multiplex` — VAE + Leiden "
        "extracts more phenotype-linked structure from the full multiplex than classical "
        "WGCNA, which dilutes the switch signal back toward the abundance baseline when "
        "abundance is added.",
        "- **Classical / multiplex WGCNA win GO-enrichment** (gene-level, "
        "abundance-biased GO). IsoGraph's low GO fraction is expected, not a failure — "
        "see the GO-invisible biology gate.",
        "- **Bottom line:** IsoGraph's defensible value is the DTU-without-DGE *content* "
        "(genes/modules structurally invisible to any abundance pipeline — biology gate, "
        "incremental association, sQTL-vs-eQTL specificity), not better module-level "
        "enrichment or phenotype rates than every WGCNA baseline. Frame the method as a "
        "complementary layer, consistent with the de-confounded gene-level result.",
    ]
    if not repl.empty:
        rsh = repl.copy()
        for c in ("mean_go_jaccard", "median_go_jaccard", "null_mean_go_jaccard"):
            if c in rsh:
                rsh[c] = rsh[c].round(4)
        if "perm_p" in rsh:
            rsh["perm_p"] = rsh["perm_p"].apply(lambda p: f"{p:.2e}" if pd.notna(p) else "NA")
        keep = [c for c in ["method", "n_preserved_aging_pairs", "n_pairs_with_go",
                            "mean_go_jaccard", "median_go_jaccard", "perm_p",
                            "null_mean_go_jaccard"] if c in rsh.columns]
        lines += [
            "",
            "## Cross-cohort GO replication consistency",
            "",
            _markdown_table(rsh[keep]),
            "",
            "Replication covers isograph vs classical `wgcna_gene` only (the matched "
            "baselines were not run through replication_go). Classical WGCNA's preserved "
            "aging modules carry more cross-cohort GO overlap — again the abundance/GO "
            "advantage — while both beat their permutation null.",
        ]
    (out_dir / "BASELINE_COMPARISON.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(description="Three-baseline IsoGraph vs WGCNA synthesis.")
    p.add_argument("--cohorts", nargs="+", default=["brainseq", "gtex"])
    p.add_argument("--fdr", type=float, default=0.10)
    args = p.parse_args()
    pooled = run(tuple(args.cohorts), args.fdr)
    if not pooled.empty:
        print(pooled.to_string(index=False))


if __name__ == "__main__":
    main()
