"""Paired statistical tests for the IsoGraph synthetic benchmark.

For each (scenario, metric) combination we compare every IsoGraph method
against WGCNA using a paired Wilcoxon signed-rank test.  Pairing is by
dataset_id / replicate seed so that each pair of measurements was produced
from identical synthetic data.

Multiple-comparison correction (Benjamini-Hochberg FDR) is applied across
all (scenario × metric × method_pair) tests jointly.

Typical usage::

    from isograph_benchmark.stats.hypothesis_tests import paired_tests
    results = paired_tests(raw_df)  # raw_df from synthetic_results.parquet
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from isograph_benchmark.stats.summarize import METRICS, benjamini_hochberg

# Primary reference method for pairwise comparisons.
REFERENCE_METHOD = "wgcna_gene"

# Columns that uniquely identify a synthetic dataset within a scenario group.
# run_scenario and run_method are excluded here because they are handled
# separately as groupby keys — including them would create duplicate columns.
_DATASET_COLS = [
    "run_dataset_id",
    "run_seed",
    "run_replicate",
    "run_n_genes",
    "run_n_samples",
    "run_switching_fraction",
    "run_noise_sd",
]


def _paired_wilcoxon(
    a: np.ndarray,
    b: np.ndarray,
) -> tuple[float, float]:
    """Return (statistic, p_value) for a paired Wilcoxon signed-rank test.

    Differences of zero are handled by scipy's default 'wilcox' zero_method.
    Returns (NaN, NaN) when fewer than 10 paired observations are available —
    the test is unreliable at very small n.
    """
    diff = a - b
    diff = diff[np.isfinite(diff)]
    if len(diff) < 10 or np.all(diff == 0):
        return np.nan, np.nan
    stat, pval = wilcoxon(diff, alternative="two-sided", zero_method="wilcox")
    return float(stat), float(pval)


def paired_tests(
    df: pd.DataFrame,
    metrics: list[str] | None = None,
    reference: str = REFERENCE_METHOD,
) -> pd.DataFrame:
    """Run paired Wilcoxon tests: each method vs. reference, per scenario × metric.

    Parameters
    ----------
    df:
        Raw results from ``synthetic_results.parquet`` (status == "completed").
    metrics:
        Metric columns to test.  Defaults to the standard METRICS list from
        ``summarize.py``.
    reference:
        Method used as the baseline for all pairwise comparisons.

    Returns
    -------
    DataFrame with columns:
        scenario, metric, method, method_ref, n_pairs,
        mean_diff (method − ref), statistic, p_value, p_adj (BH FDR),
        significant_05, significant_10, direction
    """
    if metrics is None:
        metrics = [m for m in METRICS if m in df.columns]

    scenario_col = "run_scenario" if "run_scenario" in df.columns else "scenario"
    method_col   = "run_method"   if "run_method"   in df.columns else "method"

    dataset_cols = [c for c in _DATASET_COLS if c in df.columns]

    rows: list[dict] = []

    for metric in metrics:
        if metric not in df.columns:
            continue
        sub = df[[scenario_col, method_col, metric] + dataset_cols].dropna(subset=[metric])

        for scenario, grp in sub.groupby(scenario_col):
            ref_grp = grp[grp[method_col] == reference].set_index(dataset_cols)
            if ref_grp.empty:
                continue

            methods = [m for m in grp[method_col].unique() if m != reference]
            for method in methods:
                test_grp = grp[grp[method_col] == method].set_index(dataset_cols)
                # Inner join on dataset identity — only paired observations
                common = test_grp.index.intersection(ref_grp.index)
                if len(common) == 0:
                    continue
                a = test_grp.loc[common, metric].to_numpy(float)
                b = ref_grp.loc[common, metric].to_numpy(float)
                mask = np.isfinite(a) & np.isfinite(b)
                a, b = a[mask], b[mask]

                stat, pval = _paired_wilcoxon(a, b)
                mean_diff = float(np.mean(a - b)) if len(a) > 0 else np.nan

                rows.append(
                    {
                        "scenario":   scenario,
                        "metric":     metric,
                        "method":     method,
                        "method_ref": reference,
                        "n_pairs":    int(len(a)),
                        "mean_diff":  mean_diff,
                        "statistic":  stat,
                        "p_value":    pval,
                    }
                )

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    # BH correction across all tests jointly
    result["p_adj"] = benjamini_hochberg(result["p_value"]).values
    result["significant_05"] = result["p_adj"] < 0.05
    result["significant_10"] = result["p_adj"] < 0.10
    result["direction"] = np.where(
        result["mean_diff"] > 0, "method_better",
        np.where(result["mean_diff"] < 0, "ref_better", "tied"),
    )
    return result.sort_values(["scenario", "metric", "method"]).reset_index(drop=True)
