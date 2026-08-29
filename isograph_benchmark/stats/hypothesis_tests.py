"""Paired statistical tests for the IsoGraph synthetic benchmark.

For each (scenario, metric) combination we compare every IsoGraph method
against WGCNA using a paired Wilcoxon signed-rank test.  Pairing is by
dataset_id / replicate seed so that each pair of measurements was produced
from identical synthetic data.

Effect sizes (#6).  Alongside the p-value we report two direction-aligned
effect sizes so the comparison is interpretable independent of n:

* ``rank_biserial`` — the matched-pairs rank-biserial correlation derived
  from the signed-rank statistic ((T+ − T−)/(T+ + T−)); it is the natural
  effect size for the paired Wilcoxon test, ranges [−1, 1], and is positive
  when the method beats the reference.
* ``cliffs_delta`` — Cliff's δ, a non-parametric ordinal effect size in
  [−1, 1] (P(method > ref) − P(method < ref)); reported with the conventional
  small/medium/large magnitude bins.

Multiplicity scope (#7).  All (scenario × metric × method) comparisons form a
single family.  We report Benjamini-Hochberg FDR over that *entire* family as
the canonical ``p_adj`` (the most conservative, pre-registered choice), and
additionally expose ``p_adj_within_metric`` (BH applied within each metric
across scenarios/methods) so reviewers can see how conclusions move under a
narrower family definition.  ``family`` and ``family_size`` make the
correction unit explicit on every row.

Across-scenario ranking — deliberately NOT provided.  An omnibus
Friedman/Nemenyi "critical difference" ranking over scenarios was considered
for #6 but dropped as statistically unsound here: the largest balanced block is
only ~6 scenarios (Nemenyi CD ~4 rank units on a 1-8 scale — nothing
separable), the scenarios are hand-designed stress tests rather than a random
sample of datasets, and a balanced block necessarily excludes the
confound/degradation-only methods (vae_residual, vae_reliability) that the
robustness story depends on.  The per-scenario paired tests below, reported
with rank-biserial and Cliff's-delta effect sizes, are the rigorous and
complete "how much better" comparison and replace any aggregate ranking.

Typical usage::

    from isograph_benchmark.stats.hypothesis_tests import paired_tests
    tests = paired_tests(raw_df)    # raw_df from synthetic_results.parquet
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon

from isograph_benchmark.stats.summarize import (
    METRIC_FAMILY,
    METRICS,
    benjamini_hochberg,
    bootstrap_ci,
)

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


def _signed_rank_biserial(a: np.ndarray, b: np.ndarray) -> float:
    """Matched-pairs rank-biserial correlation from the signed-rank statistic.

    r = (T+ − T−) / (T+ + T−), where T± are the sums of ranks of the positive
    and negative paired differences.  Positive => ``a`` (method) exceeds ``b``
    (reference).  Returns NaN when all differences are zero.
    """
    diff = a - b
    diff = diff[np.isfinite(diff) & (diff != 0)]
    if len(diff) == 0:
        return np.nan
    ranks = rankdata(np.abs(diff))
    t_plus = ranks[diff > 0].sum()
    t_minus = ranks[diff < 0].sum()
    total = t_plus + t_minus
    if total == 0:
        return np.nan
    return float((t_plus - t_minus) / total)


def _cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta: P(a > b) − P(a < b) over all cross-pairs, in [−1, 1]."""
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    # Vectorised pairwise sign comparison.
    diff = a[:, None] - b[None, :]
    greater = np.sum(diff > 0)
    less = np.sum(diff < 0)
    return float((greater - less) / (len(a) * len(b)))


def _cliffs_magnitude(delta: float) -> str:
    """Conventional small/medium/large bins for |Cliff's delta| (Romano 2006)."""
    if not np.isfinite(delta):
        return "na"
    d = abs(delta)
    if d < 0.147:
        return "negligible"
    if d < 0.33:
        return "small"
    if d < 0.474:
        return "medium"
    return "large"


def _paired_diff_ci(
    a: np.ndarray, b: np.ndarray, n_iter: int, alpha: float, seed: int
) -> tuple[float, float]:
    """Percentile bootstrap CI for the mean paired difference (method − reference).

    Resamples the per-simulation *differences*, which is the paired quantity — resampling
    the two arms independently would discard the pairing that makes the comparison sensitive.
    """
    diff = a - b
    diff = diff[np.isfinite(diff)]
    if len(diff) < 2:
        return np.nan, np.nan
    return bootstrap_ci(diff, n_iter=n_iter, alpha=alpha, seed=seed)


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
    n_boot: int = 10_000,
    alpha: float = 0.05,
    seed: int = 13,
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

    n_boot, alpha, seed:
        Percentile-bootstrap settings for the paired-difference confidence interval.

    Returns
    -------
    DataFrame with columns:
        scenario, metric, metric_family, method, method_ref, n_pairs,
        mean_diff (method − ref), diff_ci_low, diff_ci_high, median_diff,
        statistic, p_value, rank_biserial, cliffs_delta, cliffs_magnitude,
        family, family_size, p_adj (BH within metric_family),
        p_adj_all_metrics (BH over every test reported here),
        p_adj_within_metric (BH within each metric), significant_05,
        significant_10, direction
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
                median_diff = float(np.median(a - b)) if len(a) > 0 else np.nan
                ci_low, ci_high = _paired_diff_ci(a, b, n_iter=n_boot, alpha=alpha, seed=seed)
                rank_biserial = _signed_rank_biserial(a, b)
                cliffs = _cliffs_delta(a, b)

                rows.append(
                    {
                        "scenario":         scenario,
                        "metric":           metric,
                        "metric_family":    METRIC_FAMILY.get(metric, "core"),
                        "method":           method,
                        "method_ref":       reference,
                        "n_pairs":          int(len(a)),
                        "mean_diff":        mean_diff,
                        "median_diff":      median_diff,
                        "diff_ci_low":      ci_low,
                        "diff_ci_high":     ci_high,
                        "statistic":        stat,
                        "p_value":          pval,
                        "rank_biserial":    rank_biserial,
                        "cliffs_delta":     cliffs,
                        "cliffs_magnitude": _cliffs_magnitude(cliffs),
                    }
                )

    result = pd.DataFrame(rows)
    if result.empty:
        return result

    # ---- Multiplicity scope (#7) -------------------------------------------
    # Canonical correction: BH-FDR over every scenario × metric × method test within a
    # metric family.  The fragmentation-sensitive partition metrics are a separately
    # pre-registered addition, so they form their own family; without that split, adding 810
    # partition tests would shift every already-published `p_adj` for the original metrics by
    # enlarging a shared family.  The split bounds that specific contamination — it does not
    # freeze the core `p_adj`, which still moves as core metrics and runs accumulate.
    # `p_adj_all_metrics` keeps the pooled correction on record.
    result["family"] = "scenario_metric_method_within_metric_family"
    result["family_size"] = (
        result.groupby("metric_family")["p_value"].transform(lambda s: int(s.notna().sum()))
    )
    result["p_adj"] = (
        result.groupby("metric_family")["p_value"]
        .transform(lambda s: benjamini_hochberg(s).values)
    )
    result["p_adj_all_metrics"] = benjamini_hochberg(result["p_value"]).values

    # Sensitivity column: BH applied within each metric (a narrower family),
    # so reviewers can gauge how robust the calls are to the family definition.
    result["p_adj_within_metric"] = (
        result.groupby("metric")["p_value"]
        .transform(lambda s: benjamini_hochberg(s).values)
    )

    result["significant_05"] = result["p_adj"] < 0.05
    result["significant_10"] = result["p_adj"] < 0.10
    result["direction"] = np.where(
        result["mean_diff"] > 0, "method_better",
        np.where(result["mean_diff"] < 0, "ref_better", "tied"),
    )
    return result.sort_values(["scenario", "metric", "method"]).reset_index(drop=True)
