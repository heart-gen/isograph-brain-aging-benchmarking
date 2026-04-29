from __future__ import annotations

import numpy as np
import pandas as pd

from isograph_benchmark.config import load_yaml
from isograph_benchmark.paths import ensure_dir, rel


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


def summarize_metrics(input_path: str, metric: str = "metrics_module_recovery") -> pd.DataFrame:
    cfg = load_yaml("configs/synthetic_grid.yaml")
    df = pd.read_parquet(rel(input_path))
    n_iter = int(cfg["statistics"]["bootstrap_iterations"])
    alpha = 1 - float(cfg["statistics"]["confidence_level"])
    rows = []
    scenario_col = "run_scenario" if "run_scenario" in df.columns else "scenario"
    method_col = "run_method" if "run_method" in df.columns else "method"
    for keys, group in df.groupby([scenario_col, method_col], dropna=False):
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
    metrics = rel("benchmark", "01_synthetic", "_m", "synthetic_results.parquet")
    if not metrics.exists():
        raise SystemExit(f"Missing metrics file: {metrics}")
    out = rel("reports", "synthetic_metric_summary.parquet")
    ensure_dir(out.parent)
    summarize_metrics("benchmark/01_synthetic/_m/synthetic_results.parquet").to_parquet(
        out, index=False, compression="zstd"
    )
    print(out)


if __name__ == "__main__":
    main()
