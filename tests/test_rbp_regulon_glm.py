"""Estimability gate on the covariate-adjusted regulon GLM (reviewer item 6c).

Only module x RBP cells whose ``in_module`` coefficient has an MLE may enter the BH family.
Cells with an empty (in_module x switched) 2x2 cell, or whose fit is quasi-separated, must be
recorded with a status but no p-value, so ``multipletests`` never sees them.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data.rbp_regulon import _GLM_COVARIATES, _regulon_glm

_RNG = np.random.default_rng(13)


def _calls(switched: dict[str, list[int]], rbp: str = "RBPX") -> pd.DataFrame:
    """One region, one RBP; ``switched`` maps module_id -> per-gene 0/1 outcomes."""
    rows = []
    for module, flags in switched.items():
        for i, flag in enumerate(flags):
            rows.append({"region": "r1", "gene": f"{module}_g{i}", "module_id": module,
                         "go_invisible": False, "rbp": rbp, "switched": int(flag),
                         "pool_source": "switch_genes"})
    return pd.DataFrame(rows)


def _covariates(calls: pd.DataFrame) -> pd.DataFrame:
    genes = calls["gene"].unique()
    out = pd.DataFrame({"gene": genes})
    for col in _GLM_COVARIATES:
        out[col] = _RNG.normal(size=len(genes))
    return out


def _run(calls: pd.DataFrame) -> pd.DataFrame:
    return _regulon_glm(calls, _covariates(calls), "rbp").set_index("module_id")


def test_zero_module_cell_is_not_tested():
    """No module gene switched: ``in_module`` predicts the outcome perfectly."""
    res = _run(_calls({"M0": [0] * 20, "M1": [1] * 10 + [0] * 10,
                       "M2": [1] * 10 + [0] * 10}))
    assert res.loc["M0", "model_status"] == "separated_zero_cell"
    assert pd.isna(res.loc["M0", "pvalue_glm"])
    assert res.loc["M0", "module_switched"] == 0


def test_module_holds_every_switched_gene_is_not_tested():
    """The complement cell is empty, so the coefficient runs to +inf instead."""
    res = _run(_calls({"M0": [1] * 6 + [0] * 6, "M1": [0] * 20, "M2": [0] * 20}))
    assert res.loc["M0", "model_status"] == "separated_zero_cell"
    assert pd.isna(res.loc["M0", "pvalue_glm"])


def test_all_module_genes_switched_is_not_tested():
    res = _run(_calls({"M0": [1] * 8, "M1": [1] * 5 + [0] * 15}))
    assert res.loc["M0", "model_status"] == "separated_zero_cell"


def test_estimable_cell_is_tested():
    """Every 2x2 cell populated and plenty of events: a real, usable fit."""
    res = _run(_calls({"M0": [1] * 30 + [0] * 30, "M1": [1] * 20 + [0] * 40}))
    assert res.loc["M0", "model_status"] == "fit"
    assert np.isfinite(res.loc["M0", "pvalue_glm"])
    assert 0.0 <= res.loc["M0", "pvalue_glm"] <= 1.0
    assert np.isfinite(res.loc["M0", "odds_ratio"])


def test_bh_family_contains_only_estimable_cells():
    calls = pd.concat([
        _calls({"M0": [1] * 30 + [0] * 30, "M1": [1] * 20 + [0] * 40}, rbp="RBPA"),
        _calls({"M0": [0] * 20, "M1": [1] * 10 + [0] * 10,
                "M2": [1] * 10 + [0] * 10}, rbp="RBPB"),
    ], ignore_index=True)
    res = _regulon_glm(calls, _covariates(calls), "rbp")
    tested = res["pvalue_glm"].notna()
    assert (res.loc[tested, "model_status"] == "fit").all()
    assert res.loc[~tested, "qvalue_glm"].isna().all()
    # BH over the estimable subset only, so q >= p for every tested cell.
    assert (res.loc[tested, "qvalue_glm"] >= res.loc[tested, "pvalue_glm"] - 1e-12).all()


def test_degenerate_estimates_never_carry_a_pvalue():
    """Whatever the mechanism, a reported p-value implies a finite, plausible estimate."""
    calls = _calls({f"M{m}": list(_RNG.integers(0, 2, size=25)) for m in range(4)})
    res = _regulon_glm(calls, _covariates(calls), "rbp")
    tested = res[res["pvalue_glm"].notna()]
    assert (tested["standard_error"] < 100.0).all()
    assert (tested["coefficient"].abs() < 30.0).all()


def test_small_modules_are_skipped_before_fitting():
    res = _run(_calls({"M0": [1, 0], "M1": [1] * 20 + [0] * 20}))
    assert res.loc["M0", "model_status"] == "module_too_small"
    assert pd.isna(res.loc["M0", "pvalue_glm"])
