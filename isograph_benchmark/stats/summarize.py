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
]


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
    "run_count_dispersion",
    "run_interaction_fraction",
    "run_interaction_strength",
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

    out = rel("benchmark", "02_metrics", "_m", "synthetic_metric_summary.parquet")
    ensure_dir(out.parent)
    summary.to_parquet(out, index=False, compression="zstd")
    print(f"Wrote {len(summary):,} rows to {out.name}")

    long = make_long_metrics(df)
    out_long = rel("benchmark", "02_metrics", "_m", "synthetic_metric_long.parquet")
    long.to_parquet(out_long, index=False, compression="zstd")
    print(f"Wrote {len(long):,} rows to {out_long.name}")

    # Paired statistical tests (Wilcoxon + BH FDR) — all methods vs. WGCNA
    from isograph_benchmark.stats.hypothesis_tests import paired_tests
    print("Running paired Wilcoxon tests (vs. wgcna_gene) ...")
    tests = paired_tests(df)
    if not tests.empty:
        out_tests = rel("benchmark", "02_metrics", "_m", "synthetic_pairwise_tests.parquet")
        tests.to_parquet(out_tests, index=False, compression="zstd")
        n_sig05 = tests["significant_05"].sum()
        n_sig10 = tests["significant_10"].sum()
        print(
            f"  {len(tests):,} tests | FDR<0.05: {n_sig05} | FDR<0.10: {n_sig10}"
        )
        print(f"  Wrote {out_tests.name}")
    else:
        print("  No paired tests produced (check that wgcna_gene runs are present).")


if __name__ == "__main__":
    main()
