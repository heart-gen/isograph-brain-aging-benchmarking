from __future__ import annotations

import re

import numpy as np
import pandas as pd

from isograph_benchmark.paths import ensure_dir, rel


def _load_stats_config(yaml_path: str) -> dict:
    """Extract bootstrap_iterations and confidence_level from YAML without PyYAML."""
    from isograph_benchmark.paths import rel as _rel
    text = _rel(yaml_path).read_text()
    n_iter = int(re.search(r"bootstrap_iterations\s*:\s*(\d+)", text).group(1))
    conf = float(re.search(r"confidence_level\s*:\s*([\d.]+)", text).group(1))
    return {"bootstrap_iterations": n_iter, "confidence_level": conf}


METRICS = [
    "metrics_module_recovery",
    "metrics_switch_gene_detection_rate",
    "metrics_nonswitch_gene_module_rate",
    "metrics_n_predicted_modules",
    "metrics_n_edges",
    "measurement_elapsed_sec",
    # Abundance / isoform-role metrics (populated for the abundance scenarios;
    # missing columns are skipped automatically until a benchmark re-run emits
    # them — see compute_metrics() in benchmark/run_one.py).
    "metrics_abundance_gene_detection_rate",
    "metrics_role_switch_recall",
    "metrics_role_abundance_recall",
    "metrics_role_switch_only_n",
    "metrics_role_abundance_only_n",
    "metrics_role_coupled_n",
    "metrics_role_discordant_n",
    # A1 genetic-anchoring metrics (populated only for the genetic_anchoring scenario;
    # missing columns are skipped automatically for all other scenarios).
    "metrics_genetic_anchor_recall",
    "metrics_genetic_anchor_fpr",
    "metrics_genetic_anchor_best_r2_mean",
]

# Fragmentation-sensitive partition metrics (backfilled onto completed runs by
# benchmark.backfill_metrics; emitted natively by run_one.compute_metrics).  Best-match
# Jaccard does not penalise over-partitioning, so these are the metrics that decide whether
# a recovery advantage survives a fragmentation-aware comparison.
PARTITION_METRICS = [
    "metrics_ari_planted",                 # primary: whole-partition agreement
    "metrics_ami_planted",                 # robustness
    "metrics_ari_planted_blob",            # unassigned-gene convention sensitivity
    "metrics_v_measure_planted",
    "metrics_homogeneity_planted",         # high homogeneity + low completeness
    "metrics_completeness_planted",        #   = fragmentation; the converse = merging
    "metrics_ari_assigned",                # coverage-conditional: planted genes actually
    "metrics_ami_assigned",                #   assigned, so dropped genes are not charged
    "metrics_ari_universe",                # background genes as an extra planted class
    "metrics_ami_universe",
    "metrics_module_recovery_excess",      # observed − size-preserving null
    "metrics_module_recovery_z",
    "metrics_module_recovery_null_mean",
    "metrics_frac_planted_assigned",
]

METRICS = METRICS + PARTITION_METRICS

# Multiplicity families for hypothesis_tests.paired_tests.  The partition metrics are a
# separately pre-registered addition, so they carry their own BH family: 810 partition tests
# do not enlarge the core family and therefore cannot move a core `p_adj`.  This is not the
# same as the core `p_adj` being frozen — the core family still grows when core metrics or
# runs are added (the genetic-anchoring metrics took it from 477 to 527 tests, moving core
# `p_adj` by up to 4.8e-3), and that shift is expected.  What was checked is the part that
# matters: no FDR<0.05 call flips across the 477 previously published comparisons.
_CORE_METRICS = [m for m in METRICS if m not in PARTITION_METRICS]
METRIC_FAMILY = {
    **{m: "core" for m in _CORE_METRICS},
    **{m: "partition" for m in PARTITION_METRICS},
}


def bootstrap_ci(values: np.ndarray, n_iter: int, alpha: float, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    means = np.empty(n_iter)
    for i in range(n_iter):
        means[i] = rng.choice(values, size=len(values), replace=True).mean()
    return tuple(np.quantile(means, [alpha / 2, 1 - alpha / 2]))


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    p = pvalues.astype(float)
    valid = p.notna()
    ranked = p[valid].rank(method="first").astype(int)
    m = valid.sum()
    q = p.copy() * np.nan
    q.loc[valid] = (p[valid] * m / ranked).clip(upper=1.0)
    q.loc[valid] = q.loc[valid].sort_values(ascending=False).cummin().sort_index()
    return q


_RUN_ID_COLS = [
    "run_run_id",
    "run_replicate",
    "run_seed",
    "run_n_genes",
    "run_n_samples",
    "run_noise_sd",
    "run_switching_fraction",
    "run_abundance_imbalance",
    "run_abundance_fraction",
    "run_count_dispersion",
    "run_interaction_fraction",
    "run_interaction_strength",
    # Gap #6 confound-scenario sweep parameters + interpretation/degradation knobs;
    # retained so the long table can facet the confound-robustness, degradation
    # fallback, and multi-isoform interpretation figures.
    "run_degradation_3p_bias",
    "run_cell_composition_cv",
    "run_batch_effect_sd",
    "run_library_depth_cv",
    "run_dual_signal_fraction",
    "run_n_transcripts_per_gene",
]


def make_long_metrics(df: pd.DataFrame) -> pd.DataFrame:
    scenario_col = "run_scenario" if "run_scenario" in df.columns else "scenario"
    method_col = "run_method" if "run_method" in df.columns else "method"
    id_cols = [scenario_col, method_col] + [c for c in _RUN_ID_COLS if c in df.columns]
    available = [m for m in METRICS if m in df.columns]
    long = (
        df[id_cols + available]
        .melt(id_vars=id_cols, value_vars=available, var_name="metric", value_name="value")
        .rename(columns={scenario_col: "scenario", method_col: "method"})
        .reset_index(drop=True)
    )
    return long


def summarize_metrics(
    df: pd.DataFrame,
    metric: str,
    n_iter: int,
    alpha: float,
) -> pd.DataFrame:
    scenario_col = "run_scenario" if "run_scenario" in df.columns else "scenario"
    method_col = "run_method" if "run_method" in df.columns else "method"
    rows = []
    for keys, group in df.groupby([scenario_col, method_col], dropna=False):
        if metric not in group.columns:
            continue
        vals = group[metric].to_numpy(dtype=float)
        lo, hi = bootstrap_ci(vals, n_iter=n_iter, alpha=alpha)
        rows.append(
            {
                "scenario": keys[0],
                "method": keys[1],
                "metric": metric,
                "n": int(np.isfinite(vals).sum()),
                "mean": float(np.nanmean(vals)),
                "median": float(np.nanmedian(vals)),
                "ci_low": lo,
                "ci_high": hi,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    results_path = rel("benchmark", "01_synthetic", "_m", "synthetic_results.parquet")
    if not results_path.exists():
        raise SystemExit(f"Missing results file: {results_path}\nRun collect_results first.")

    cfg = _load_stats_config("configs/synthetic_grid.yaml")
    n_iter = int(cfg["bootstrap_iterations"])
    alpha = 1 - float(cfg["confidence_level"])

    df = pd.read_parquet(results_path)
    df = df[df["status"] == "completed"].copy()
    print(f"Summarizing {len(df):,} completed runs across {len(METRICS)} metrics")

    parts = []
    for metric in METRICS:
        if metric not in df.columns:
            print(f"  Skipping {metric} (column not found)")
            continue
        part = summarize_metrics(df, metric=metric, n_iter=n_iter, alpha=alpha)
        parts.append(part)
        print(f"  {metric}: {len(part)} group combinations")

    summary = pd.concat(parts, ignore_index=True)

    out = rel("benchmark", "03_metrics", "_m", "synthetic_metric_summary.parquet")
    ensure_dir(out.parent)
    summary.to_parquet(out, index=False, compression="zstd")
    print(f"Wrote {len(summary):,} rows to {out.name}")

    long = make_long_metrics(df)
    out_long = rel("benchmark", "03_metrics", "_m", "synthetic_metric_long.parquet")
    long.to_parquet(out_long, index=False, compression="zstd")
    print(f"Wrote {len(long):,} rows to {out_long.name}")

    # Paired statistical tests (Wilcoxon + effect sizes + BH FDR) vs. WGCNA
    from isograph_benchmark.stats.hypothesis_tests import paired_tests
    print("Running paired Wilcoxon tests + paired-difference bootstrap CIs (vs. wgcna_gene) ...")
    tests = paired_tests(df, n_boot=n_iter, alpha=alpha, seed=13)
    if not tests.empty:
        out_tests = rel("benchmark", "03_metrics", "_m", "synthetic_pairwise_tests.parquet")
        tests.to_parquet(out_tests, index=False, compression="zstd")
        n_sig05 = tests["significant_05"].sum()
        n_sig10 = tests["significant_10"].sum()
        n_large = (tests["cliffs_magnitude"] == "large").sum()
        fams = tests.groupby("metric_family")["family_size"].first().to_dict()
        print(
            f"  {len(tests):,} tests (BH families={fams}) | FDR<0.05: {n_sig05} | "
            f"FDR<0.10: {n_sig10} | large |Cliff's d|: {n_large}"
        )
        print(f"  Wrote {out_tests.name}")
    else:
        print("  No paired tests produced (check that wgcna_gene runs are present).")


if __name__ == "__main__":
    main()
