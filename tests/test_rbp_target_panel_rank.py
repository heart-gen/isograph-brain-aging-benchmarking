"""Ranking modes for the RBP perturbation panel (reviewer item 6c, task #12).

`--rank-by adjusted` orders candidates on the covariate-adjusted regulon GLM instead of the
raw hypergeometric enrichment. The ordering key is the *lower bound* of the adjusted 95% CI,
so a large but imprecise odds ratio cannot outrank a smaller one whose interval excludes 1,
and a module the adjustment contradicts sinks below every module it supports.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data import rbp_target_panel as rtp
from isograph_benchmark.real_data.rbp_target_panel import _best_module, rank_candidates


@pytest.fixture(autouse=True)
def _no_coloc(tmp_path, monkeypatch):
    """Keep the fixtures hermetic: no deep-dive parquet, so no coloc anchors are joined."""
    monkeypatch.setattr(rtp, "_DEEP_DIVE", tmp_path / "absent")


def _ev(rows: list[dict]) -> pd.DataFrame:
    """Build a (gene, region) evidence frame; `rows` overrides the neutral defaults."""
    base = {"region": "r1", "module_id": "M0", "module_go": "go", "go_invisible": False,
            "pool_source": "switch_genes", "pheno_module": True, "binding_supported": True,
            "switch_active": True, "switch_r": 0.5, "module_role": "driver",
            "pheno_fdr": 0.01, "module_size": 100, "n_switched": 50,
            "module_q": 1e-3, "module_enrichment": 1.5, "module_q_glm": 0.5,
            "module_or_glm": 1.0, "module_or_ci_low": 0.8, "module_or_ci_high": 1.3}
    return pd.DataFrame([{**base, **r} for r in rows])


def _order(cand: pd.DataFrame) -> list[str]:
    return list(cand.sort_values("rank")["gene"])


def test_adjusted_prefers_the_tighter_interval_over_the_bigger_point_estimate():
    """OR 1.60 (1.09-2.35) must outrank OR 1.85 (0.96-3.53): the second interval covers 1."""
    ev = _ev([
        {"gene": "BIG", "module_id": "M1", "module_or_glm": 1.85,
         "module_or_ci_low": 0.96, "module_or_ci_high": 3.53},
        {"gene": "TIGHT", "module_id": "M2", "module_or_glm": 1.60,
         "module_or_ci_low": 1.09, "module_or_ci_high": 2.35},
    ])
    assert _order(rank_candidates(ev, rank_by="adjusted")) == ["TIGHT", "BIG"]


def test_adjusted_sinks_a_module_the_adjustment_contradicts():
    """Depleted after adjustment (CI wholly below 1) loses to a supported module, even with
    stronger unadjusted enrichment and wider regional recurrence."""
    ev = _ev([
        {"gene": "DEPLETED", "region": "r1", "module_id": "M1", "module_q": 1e-8,
         "module_enrichment": 1.42, "module_or_glm": 0.56,
         "module_or_ci_low": 0.38, "module_or_ci_high": 0.82},
        {"gene": "DEPLETED", "region": "r2", "module_id": "M1", "module_q": 1e-8,
         "module_enrichment": 1.42, "module_or_glm": 0.56,
         "module_or_ci_low": 0.38, "module_or_ci_high": 0.82},
        {"gene": "SUPPORTED", "region": "r1", "module_id": "M2", "module_q": 1e-2,
         "module_enrichment": 1.10, "module_or_glm": 1.60,
         "module_or_ci_low": 1.09, "module_or_ci_high": 2.35},
    ])
    cand = rank_candidates(ev, rank_by="adjusted").set_index("gene")
    assert cand.loc["SUPPORTED", "rank"] < cand.loc["DEPLETED", "rank"]
    assert bool(cand.loc["SUPPORTED", "adjusted_or_ok"])
    assert not bool(cand.loc["DEPLETED", "adjusted_or_ok"])
    # ... and the reverse under the pre-adjustment rule, which sees only recurrence.
    old = rank_candidates(ev, rank_by="hypergeometric").set_index("gene")
    assert old.loc["DEPLETED", "rank"] < old.loc["SUPPORTED", "rank"]


def test_hypergeometric_mode_ignores_the_adjusted_columns():
    """The legacy ordering must not shift when only the adjusted arm changes."""
    rows = [
        {"gene": "A", "region": "r1", "module_id": "M1", "module_q": 1e-9,
         "module_enrichment": 2.0},
        {"gene": "B", "region": "r1", "module_id": "M2", "module_q": 1e-2,
         "module_enrichment": 1.1},
    ]
    plain = rank_candidates(_ev(rows), rank_by="hypergeometric")
    flipped = [{**r, "module_or_glm": 9.0, "module_or_ci_low": 5.0, "module_or_ci_high": 15.0}
               if r["gene"] == "B" else r for r in rows]
    assert _order(plain) == _order(rank_candidates(_ev(flipped), rank_by="hypergeometric"))
    # under the adjusted rule the same change does move B ahead
    assert _order(rank_candidates(_ev(flipped), rank_by="adjusted")) == ["B", "A"]


def test_min_adjusted_or_threshold_is_applied():
    ev = _ev([
        {"gene": "WEAK", "module_id": "M1", "module_or_glm": 1.2,
         "module_or_ci_low": 1.05, "module_or_ci_high": 1.4},
        {"gene": "STRONG", "module_id": "M2", "module_or_glm": 2.5,
         "module_or_ci_low": 1.10, "module_or_ci_high": 5.0},
    ])
    permissive = rank_candidates(ev, rank_by="adjusted", min_adjusted_or=1.0)
    assert set(permissive["adjusted_or_ok"]) == {True}
    strict = rank_candidates(ev, rank_by="adjusted", min_adjusted_or=2.0).set_index("gene")
    assert bool(strict.loc["STRONG", "adjusted_or_ok"])
    assert not bool(strict.loc["WEAK", "adjusted_or_ok"])
    assert strict.loc["STRONG", "rank"] < strict.loc["WEAK", "rank"]


def test_best_module_differs_between_the_two_modes():
    """A gene in two modules is reported under the one its ranking mode considers strongest."""
    ev = _ev([
        {"gene": "G", "region": "rHYP", "module_id": "M_hyper", "module_q": 1e-20,
         "module_enrichment": 2.0, "module_or_glm": 0.95,
         "module_or_ci_low": 0.78, "module_or_ci_high": 1.18},
        {"gene": "G", "region": "rADJ", "module_id": "M_adj", "module_q": 1e-2,
         "module_enrichment": 1.2, "module_or_glm": 1.60,
         "module_or_ci_low": 1.09, "module_or_ci_high": 2.35},
    ])
    assert _best_module(ev, "hypergeometric")["best_module"].iloc[0] == "M_hyper"
    assert _best_module(ev, "adjusted")["best_module"].iloc[0] == "M_adj"


def test_adjusted_ci_flag_marks_only_intervals_excluding_one():
    ev = _ev([
        {"gene": "EXCLUDES", "module_id": "M1", "module_or_glm": 1.6,
         "module_or_ci_low": 1.09, "module_or_ci_high": 2.35},
        {"gene": "COVERS", "module_id": "M2", "module_or_glm": 1.85,
         "module_or_ci_low": 0.96, "module_or_ci_high": 3.53},
    ])
    cand = rank_candidates(ev, rank_by="adjusted").set_index("gene")
    assert bool(cand.loc["EXCLUDES", "adjusted_ci_excludes_null"])
    assert not bool(cand.loc["COVERS", "adjusted_ci_excludes_null"])


def test_ineligible_genes_sort_last_in_both_modes():
    ev = _ev([
        {"gene": "ELIGIBLE", "module_id": "M1", "module_or_glm": 1.1,
         "module_or_ci_low": 1.01, "module_or_ci_high": 1.3},
        {"gene": "NO_PHENO", "module_id": "M2", "pheno_module": False,
         "pool_source": "all_genes", "module_q": 1e-30, "module_enrichment": 5.0,
         "module_or_glm": 9.0, "module_or_ci_low": 5.0, "module_or_ci_high": 15.0},
    ])
    for mode in ("adjusted", "hypergeometric"):
        assert _order(rank_candidates(ev, rank_by=mode)) == ["ELIGIBLE", "NO_PHENO"], mode


def test_missing_adjusted_columns_fall_back_without_error():
    """Older stage-2 tables have no GLM arm; the CLI must still rank."""
    ev = _ev([{"gene": "A", "module_id": "M1"}, {"gene": "B", "module_id": "M2"}]).drop(
        columns=["module_or_glm", "module_or_ci_low", "module_or_ci_high", "module_q_glm"])
    cand = rank_candidates(ev, rank_by="adjusted")
    assert len(cand) == 2
    assert cand["best_module_q_glm"].isna().all()


def test_ranking_is_deterministic():
    ev = _ev([{"gene": g, "module_id": "M1"} for g in ("C", "A", "B")])
    first = _order(rank_candidates(ev, rank_by="adjusted"))
    assert first == _order(rank_candidates(ev.sample(frac=1, random_state=13),
                                           rank_by="adjusted"))
    assert first == ["A", "B", "C"]   # symbol is the final tie-break


def test_best_module_estimates_come_from_one_module():
    """The reported OR and its CI must be the selected module's own triple.

    Aggregating each column independently across the gene's modules (a separate max per
    column) produced intervals no model fitted -- e.g. a point estimate and upper bound from
    one module against a lower bound from another.
    """
    ev = _ev([
        # the module the adjusted rule selects: highest lower bound, modest OR, tight top
        {"gene": "G", "region": "rSEL", "module_id": "M_sel", "module_enrichment": 1.1,
         "module_q": 1e-2, "module_or_glm": 1.55, "module_or_ci_low": 1.03,
         "module_or_ci_high": 2.60},
        # a rival carrying the larger OR and the larger upper bound, but a lower bound under 1
        {"gene": "G", "region": "rBIG", "module_id": "M_big", "module_enrichment": 2.0,
         "module_q": 1e-9, "module_or_glm": 1.85, "module_or_ci_low": 0.96,
         "module_or_ci_high": 3.53},
    ])
    r = rank_candidates(ev, rank_by="adjusted").iloc[0]
    assert r["best_module"] == "M_sel"
    assert r["best_module_or"] == pytest.approx(1.55)
    assert r["best_module_or_ci_low"] == pytest.approx(1.03)
    assert r["best_module_or_ci_high"] == pytest.approx(2.60)   # not 3.53, from M_big
    # and the interval must contain its own point estimate
    assert r["best_module_or_ci_low"] <= r["best_module_or"] <= r["best_module_or_ci_high"]
