"""Driver selection and distinct-pair pooling for the long-read driver confirmation."""
import numpy as np
import pandas as pd

from isograph_benchmark.real_data.longread_driver_confirm import (
    distinct_pairs,
    flag_driver_pairs,
    top_k_drivers,
)


def _pairs():
    cols = ["region", "gene", "t1", "t2", "t1_detected", "t2_detected", "pair_detected",
            "t1_mean_if", "t2_mean_if", "usage_spearman", "switch_like_strict"]
    rows = [
        # same pair nominated in both regions, identical long-read outcome
        ("frontal_cortex_ba9", "G1", "T1", "T2", True, True, True, 0.40, 0.30, -0.5, True),
        ("cortex", "G1", "T2", "T1", True, True, True, 0.30, 0.40, -0.5, True),
        # driver member (T3) at 1% of gene output: not usable; T6 is not a driver
        ("cortex", "G1", "T6", "T3", True, True, True, 0.40, 0.01, -0.2, False),
        # no driver member
        ("cortex", "G2", "T4", "T5", True, False, False, 0.20, 0.00, np.nan, False),
    ]
    p = pd.DataFrame(rows, columns=cols)
    p["max_mean_if"] = p[["t1_mean_if", "t2_mean_if"]].max(axis=1)
    p["switch_like"] = p["pair_detected"] & (p["usage_spearman"] < 0)
    return p


def _drivers():
    return pd.DataFrame({
        "region": ["frontal_cortex_ba9", "cortex", "cortex"],
        "module_id": ["M000", "M003", "M004"],
        "transcript_id": ["T1", "T1", "T3"],
    })


def test_top_k_drivers_ranks_by_abs_r_and_strips_versions():
    pol = pd.DataFrame({"transcript_id": ["ENST1.2", "ENST2.1", "ENST3.5", "ENST4.1"],
                        "gene_id": ["ENSG1.3"] * 4, "r": [0.1, -0.9, 0.5, np.nan]})
    d = top_k_drivers(pol, 2)
    assert list(d["transcript_id"]) == ["ENST2", "ENST3"]
    assert list(d["rank"]) == [1, 2] and d["gene"].iloc[0] == "ENSG1"


def test_flag_driver_pairs_uses_the_driver_side_abundance():
    f = flag_driver_pairs(_pairs(), _drivers(), min_if=0.05)
    assert f["driver_pair"].tolist() == [True, True, True, False]
    # T3 is the driver in row 2 and sits at 1% -> not usable; T1 rows are usable
    assert f["driver_usable"].tolist() == [True, True, False, False]
    assert f.loc[2, "driver_mean_if"] == 0.01
    assert f.loc[1, "driver_modules"] == "M003"


def test_distinct_pairs_counts_a_pair_once_across_regions():
    dp = distinct_pairs(flag_driver_pairs(_pairs(), _drivers(), min_if=0.05))
    assert len(dp) == 3
    both = dp[(dp["tx_a"] == "T1") & (dp["tx_b"] == "T2")].iloc[0]
    assert both["n_regions"] == 2 and both["regions"] == "cortex;frontal_cortex_ba9"
    assert both["driver_modules"] == "M000;M003" and both["switch_like"]
    assert dp[dp["gene"] == "G2"]["driver_pair"].iloc[0] is np.False_
