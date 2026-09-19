"""Pairwise `direction` must read each metric against its polarity."""
import numpy as np
import pandas as pd

from isograph_benchmark.stats.hypothesis_tests import paired_tests
from isograph_benchmark.stats.summarize import METRIC_POLARITY, METRICS


def _frame(metric: str, method_value: float, ref_value: float, n: int = 12) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    rows = []
    for i in range(n):
        noise = rng.normal(0, 0.01)
        for method, value in (("isograph_vae", method_value), ("wgcna_gene", ref_value)):
            rows.append({"run_scenario": "s", "run_method": method, "run_dataset_id": f"d{i}",
                         "run_seed": i, metric: value + noise + (0.001 * i if method == "isograph_vae" else 0)})
    return pd.DataFrame(rows)


def test_every_metric_has_a_polarity():
    assert set(METRICS) <= set(METRIC_POLARITY)
    assert set(METRIC_POLARITY.values()) <= {-1, 0, 1}


def test_lower_is_better_metric_reads_a_lower_value_as_better():
    metric = "metrics_nonswitch_gene_module_rate"
    out = paired_tests(_frame(metric, 0.1, 0.5), metrics=[metric], n_boot=50)
    row = out.iloc[0]
    assert row["mean_diff"] < 0
    assert row["direction_raw"] == "method_lower"
    assert row["direction"] == "method_better"
    assert row["higher_is_better"] is False or row["higher_is_better"] == False  # noqa: E712


def test_higher_is_better_metric_keeps_the_raw_sign():
    metric = "metrics_module_recovery"
    out = paired_tests(_frame(metric, 0.9, 0.5), metrics=[metric], n_boot=50)
    assert out.iloc[0]["direction"] == "method_better"


def test_descriptive_metric_has_no_better_side():
    metric = "metrics_n_predicted_modules"
    out = paired_tests(_frame(metric, 40.0, 10.0), metrics=[metric], n_boot=50)
    row = out.iloc[0]
    assert row["direction"] == "not_applicable"
    assert pd.isna(row["higher_is_better"])
