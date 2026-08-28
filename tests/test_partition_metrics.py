"""Unit tests for the fragmentation-sensitive partition metrics (reviewer items 1 and 4)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph.evaluation.metrics import module_recovery_score

from isograph_benchmark.benchmark.partition_metrics import (
    best_match_jaccard,
    partition_labels,
    partition_metrics,
    size_preserving_null,
)


def _frame(assignments: dict[str, list[str]]) -> pd.DataFrame:
    rows = [(gene, mod) for mod, genes in assignments.items() for gene in genes]
    return pd.DataFrame(rows, columns=["gene_id", "module_id"])


@pytest.fixture
def truth() -> pd.DataFrame:
    return _frame({
        "T0": [f"g{i}" for i in range(10)],
        "T1": [f"g{i}" for i in range(10, 20)],
    })


@pytest.fixture
def universe() -> list[str]:
    return [f"g{i}" for i in range(60)]


# --------------------------------------------------------------------------- #
# Item 1: whole-partition metrics behave as advertised
# --------------------------------------------------------------------------- #
def test_perfect_recovery(truth, universe):
    m = partition_metrics(truth.copy(), truth, universe=universe, n_perm=50)
    assert m["ari_planted"] == pytest.approx(1.0)
    assert m["ami_planted"] == pytest.approx(1.0)
    assert m["v_measure_planted"] == pytest.approx(1.0)
    assert m["module_recovery_recomputed"] == pytest.approx(1.0)


def test_fragmentation_is_high_homogeneity_low_completeness(truth, universe):
    """Splitting every truth module into singletons: pure but maximally incomplete."""
    shattered = _frame({f"P{i}": [g] for i, g in enumerate(truth["gene_id"])})
    m = partition_metrics(shattered, truth, universe=universe, n_perm=50)
    assert m["homogeneity_planted"] == pytest.approx(1.0)
    assert m["completeness_planted"] < 0.6
    assert m["ari_planted"] < 0.05


def test_merging_is_the_converse(truth, universe):
    """Collapsing both truth modules into one: complete but impure."""
    merged = _frame({"P0": list(truth["gene_id"])})
    m = partition_metrics(merged, truth, universe=universe, n_perm=50)
    assert m["completeness_planted"] == pytest.approx(1.0)
    assert m["homogeneity_planted"] < 0.05


def test_random_labels_give_ari_near_zero(truth, universe):
    rng = np.random.default_rng(0)
    genes = list(truth["gene_id"])
    mods = rng.integers(0, 2, size=len(genes))
    random_pred = pd.DataFrame({"gene_id": genes, "module_id": [f"P{m}" for m in mods]})
    m = partition_metrics(random_pred, truth, universe=universe, n_perm=50)
    assert abs(m["ari_planted"]) < 0.35


def test_unassigned_conventions_differ(truth, universe):
    """Singleton vs pooled unassigned genes are different partitions, so both are reported."""
    partial = _frame({"P0": [f"g{i}" for i in range(5)], "P1": [f"g{i}" for i in range(10, 15)]})
    m = partition_metrics(partial, truth, universe=universe, n_perm=50)
    assert m["ari_planted"] != pytest.approx(m["ari_planted_blob"])
    assert m["frac_planted_assigned"] == pytest.approx(0.5)


def test_coverage_conditional_ari_ignores_dropped_genes(truth, universe):
    """Half of each truth module assigned, and assigned perfectly.

    ``ari_planted`` is dragged down by the genes the method never placed; ``ari_assigned``
    scores only what it did place, so the two together separate a coverage shortfall from a
    grouping error.
    """
    partial = _frame({"P0": [f"g{i}" for i in range(5)], "P1": [f"g{i}" for i in range(10, 15)]})
    m = partition_metrics(partial, truth, universe=universe, n_perm=50)
    assert m["ari_assigned"] == pytest.approx(1.0)
    assert m["ami_assigned"] == pytest.approx(1.0)
    assert m["ari_planted"] < 1.0


def test_coverage_conditional_ari_still_penalises_mis_grouping(truth, universe):
    """Full coverage but the two truth modules interleaved: coverage is not the problem."""
    genes = [f"g{i}" for i in range(20)]
    scrambled = _frame({"P0": genes[::2], "P1": genes[1::2]})
    m = partition_metrics(scrambled, truth, universe=universe, n_perm=50)
    assert m["frac_planted_assigned"] == pytest.approx(1.0)
    assert m["ari_assigned"] < 0.05
    assert m["ari_assigned"] == pytest.approx(m["ari_planted"])


def test_single_gene_truth_suppresses_the_partition_metrics(universe):
    """``negative_control_noise`` plants one gene in one module.

    Every method then scores a vacuous ari/homogeneity/completeness of 1.0 on it, which is
    exactly the free score these metrics exist to remove, so the comparison metrics must be
    withheld. The size-preserving null is still well-defined and must still be reported —
    it is what actually diagnoses that scenario.
    """
    tiny = _frame({"T0": ["g0"]})
    pred = _frame({f"P{i}": [f"g{i}"] for i in range(20)})
    m = partition_metrics(pred, tiny, universe=universe, n_perm=50)
    assert m["ari_planted"] is None
    assert m["homogeneity_planted"] is None
    assert m["ari_assigned"] is None
    assert m["n_truth_modules"] == pytest.approx(1.0)
    assert m["module_recovery_recomputed"] is not None
    assert m["module_recovery_perm_p"] is not None


def test_partition_labels_cover_the_whole_universe(truth, universe):
    t, p = partition_labels(truth.copy(), truth, universe=universe)
    assert len(t) == len(p) == len(universe)
    # Background genes collapse to a single truth class; planted genes keep their own.
    assert len(np.unique(t)) == 3


# --------------------------------------------------------------------------- #
# Item 4: the null must match the published metric's conventions
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_best_match_jaccard_matches_published_metric(truth, seed):
    """Our reimplementation must equal module_recovery_score on arbitrary partitions,
    including when predicted modules contain background genes absent from truth."""
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(60)]
    n_mod = int(rng.integers(2, 12))
    pred = pd.DataFrame({
        "gene_id": genes,
        "module_id": [f"P{m}" for m in rng.integers(0, n_mod, size=len(genes))],
    })
    assert best_match_jaccard(pred, truth) == pytest.approx(module_recovery_score(pred, truth))


def test_null_preserves_the_predicted_module_size_multiset(truth):
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(60)]
    pred = pd.DataFrame({
        "gene_id": genes,
        "module_id": [f"P{m}" for m in rng.integers(0, 5, size=len(genes))],
    })
    sizes = sorted(pred.groupby("module_id").size())
    # The permutation only reorders the label vector, so sizes are invariant by construction;
    # assert it explicitly because the whole calibration depends on it.
    codes = pd.factorize(pred["module_id"], sort=True)[0]
    permuted = np.random.default_rng(1).permutation(codes)
    assert sorted(np.bincount(permuted)) == sizes


def test_fragmenting_method_gets_a_large_free_recovery(truth, universe):
    """The signal-free WGCNA pathology: many tiny modules score well above zero on
    best-match Jaccard, and the size-preserving null is what exposes it."""
    rng = np.random.default_rng(7)
    many = pd.DataFrame({
        "gene_id": universe,
        "module_id": [f"P{m}" for m in rng.integers(0, 25, size=len(universe))],
    })
    m = partition_metrics(many, truth, universe=universe, n_perm=200, seed=13)
    assert m["module_recovery_null_mean"] > 0.05          # recovery is not free-of-charge zero
    assert m["module_recovery_perm_p"] > 0.05             # random labels are not significant
    assert abs(m["module_recovery_z"]) < 3


def test_null_is_deterministic_under_a_fixed_seed(truth, universe):
    rng = np.random.default_rng(3)
    pred = pd.DataFrame({
        "gene_id": universe,
        "module_id": [f"P{m}" for m in rng.integers(0, 6, size=len(universe))],
    })
    a = size_preserving_null(pred, truth, n_perm=100, seed=13)
    b = size_preserving_null(pred, truth, n_perm=100, seed=13)
    assert np.array_equal(a, b)


def test_empty_predictions_return_none_not_zero(truth, universe):
    empty = pd.DataFrame(columns=["gene_id", "module_id"])
    m = partition_metrics(empty, truth, universe=universe, n_perm=10)
    assert all(v is None for v in m.values())
