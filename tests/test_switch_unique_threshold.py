"""Threshold sensitivity: the category rule, and calls that must survive every alpha."""
import numpy as np
import pandas as pd
import pytest

sut = pytest.importorskip("isograph_benchmark.real_data.switch_unique_threshold")


def _gl(sw, ab):
    return pd.DataFrame({"gene_id": [f"g{i}" for i in range(len(sw))],
                         "fdr_switch_given_abund": sw, "fdr_abund_given_switch": ab})


def test_category_rule_matches_incremental_association_including_the_boundary():
    # Boundaries sit ON alpha: incremental_association uses `<=`, so 0.10 counts as
    # significant at alpha = 0.10. A gene at the boundary on the abundance side is
    # therefore "both", not switch-unique.
    gl = _gl(sw=[0.01, 0.10, 0.10, 0.50, 0.01],
             ab=[0.50, 0.50, 0.10, 0.50, 0.01])
    cat = sut._categorize(gl, 0.10).tolist()
    assert cat == ["switch_unique", "switch_unique", "both", "neither", "both"]

    # Tightening to 0.05 moves the 0.10 genes out of significance entirely.
    assert sut._categorize(gl, 0.05).tolist() == [
        "switch_unique", "neither", "neither", "neither", "both"]


def test_a_call_needs_every_alpha_to_agree():
    df = pd.DataFrame({
        "label": ["collapse"] * 3 + ["retain"] * 3 + ["flips"] * 3,
        "alpha": [0.05, 0.10, 0.20] * 3,
        "switch_unique_base": [100] * 9,
        "persistence": [0.0, 0.01, 0.02,       # below COLLAPSE_MAX at every alpha
                        0.30, 0.40, 0.50,      # at or above RETAIN_MIN at every alpha
                        0.01, 0.40, 0.30],     # collapses at one alpha, retains at others
    })
    calls = sut._robust_calls(df).set_index("label")["call"].to_dict()
    assert calls["collapse"] == "collapses at every alpha"
    assert calls["retain"] == "retains at every alpha"
    assert calls["flips"] == "threshold-dependent"


def test_a_call_resting_on_a_handful_of_genes_is_flagged():
    df = pd.DataFrame({
        "label": ["tiny"] * 3,
        "alpha": [0.05, 0.10, 0.20],
        "switch_unique_base": [2, 2, 5],       # below MIN_DENOMINATOR
        "persistence": [0.5, 1.0, 0.8],
    })
    call = sut._robust_calls(df)["call"].iloc[0]
    assert call.endswith("(small set)")


def test_rank_stability_compares_each_alpha_with_the_reference():
    ref = sut.REFERENCE_ALPHA
    labels = [f"a{i}" for i in range(6)]
    df = pd.concat([
        pd.DataFrame({"label": labels, "alpha": ref,
                      "switch_unique_base": [10, 20, 30, 40, 50, 60],
                      "persistence": np.linspace(0.1, 0.6, 6)}),
        # A uniform rescale preserves every rank.
        pd.DataFrame({"label": labels, "alpha": 0.20,
                      "switch_unique_base": [20, 40, 60, 80, 100, 120],
                      "persistence": np.linspace(0.2, 0.7, 6)}),
    ], ignore_index=True)
    ranks = sut._rank_stability(df).set_index("alpha")

    assert ranks.loc[ref, "spearman_switch_unique_base"] == pytest.approx(1.0)
    assert ranks.loc[0.20, "spearman_switch_unique_base"] == pytest.approx(1.0)
    assert ranks.loc[0.20, "spearman_persistence"] == pytest.approx(1.0)
