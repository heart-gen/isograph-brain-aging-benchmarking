"""Tests for the held-out module-context test of DTU evidence (PI review item 12a)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data import module_dtu_added_value as m


# --------------------------------------------------------------------------- #
# per-gene evidence
# --------------------------------------------------------------------------- #
def test_simes_matches_hand_computation():
    res = pd.DataFrame({
        "gene_id": ["G1.1", "G1.1", "G1.1", "G2.3", "G2.3"],
        "pval": [0.04, 0.01, 0.30, 0.5, np.nan],
    })
    out = m.simes_by_gene(res).set_index("gene")
    # G1: sorted 0.01, 0.04, 0.30 -> min(3*.01/1, 3*.04/2, 3*.30/3) = 0.03
    assert out.loc["G1", "simes_p"] == pytest.approx(0.03)
    assert out.loc["G1", "min_p"] == pytest.approx(0.01)
    assert out.loc["G1", "n_tx"] == 3
    # G2: the NaN isoform is dropped, not counted
    assert out.loc["G2", "n_tx"] == 1
    assert out.loc["G2", "simes_p"] == pytest.approx(0.5)
    assert out.loc["G1", "evidence"] == pytest.approx(-np.log10(0.03))


def test_gene_covariates_minor_usage():
    counts = np.array([[90, 80], [10, 20], [5, 5]], dtype=float)  # G1 two isoforms, G2 one
    out = m.gene_covariates(counts, pd.Series(["G1.1", "G1.1", "G2.1"])).set_index("gene")
    assert out.loc["G1", "minor_usage"] == pytest.approx(1 - np.mean([0.9, 0.8]))
    assert out.loc["G2", "minor_usage"] == pytest.approx(0.0)
    assert out.loc["G1", "log_mean_count"] == pytest.approx(np.log1p(100.0))


# --------------------------------------------------------------------------- #
# context and permutation machinery
# --------------------------------------------------------------------------- #
def test_loo_context_excludes_self():
    labels = np.array(["a", "a", "a", "b", "b"])
    vals = np.array([1.0, 2.0, 6.0, 10.0, 20.0])
    ctx = m.loo_context(labels, vals)
    assert ctx == pytest.approx([4.0, 3.5, 1.5, 20.0, 10.0])


def test_loo_context_matrix_matches_vector_version():
    rng = np.random.default_rng(0)
    labels = rng.integers(0, 7, size=200)
    vals = rng.normal(size=200)
    mat = np.column_stack([rng.permutation(labels) for _ in range(5)])
    got = m.loo_context_matrix(mat, vals)
    for j in range(5):
        assert got[:, j] == pytest.approx(m.loo_context(mat[:, j], vals))


def test_permute_within_preserves_each_label_stratum_makeup():
    rng = np.random.default_rng(1)
    labels = rng.integers(0, 6, size=500)
    strata = rng.integers(0, 9, size=500)
    perm = m.permute_within(labels, strata, rng)
    for s in np.unique(strata):
        assert sorted(perm[strata == s]) == sorted(labels[strata == s])
    assert not np.array_equal(perm, labels)


def test_split_indices_are_the_stage04_halves():
    from isograph_benchmark.real_data import stability
    for seed in (0, 3):
        a, b = stability._split_indices(101, stability.SEED_BASE + seed)
        assert np.array_equal(m.split_indices(101, seed, "A"), a)
        assert np.array_equal(m.split_indices(101, seed, "B"), b)
        assert set(a).isdisjoint(b) and len(a) + len(b) == 101


def test_design_matrix_is_full_rank():
    rng = np.random.default_rng(2)
    d = pd.DataFrame({"e_disc": rng.exponential(size=300), "n_tx": rng.integers(2, 9, 300),
                      "log_mean_count": rng.normal(5, 1, 300), "minor_usage": rng.uniform(0, .5, 300)})
    X = m.design_matrix(d)
    assert np.linalg.matrix_rank(X) == X.shape[1]


def test_residualizer_handles_rank_deficiency():
    rng = np.random.default_rng(3)
    x = rng.normal(size=100)
    X = np.column_stack([np.ones(100), x, 2 * x])       # rank 2
    v = rng.normal(size=100)
    r = m._residualizer(X)(v)
    expected = v - np.column_stack([np.ones(100), x]) @ np.linalg.lstsq(
        np.column_stack([np.ones(100), x]), v, rcond=None)[0]
    assert r == pytest.approx(expected)


# --------------------------------------------------------------------------- #
# the test itself, on simulated halves
# --------------------------------------------------------------------------- #
def _sim(n_genes=3000, n_modules=30, module_signal=1.0, technical=0.0, seed=0):
    """Two halves of per-gene evidence around a shared true signal theta.

    ``module_signal``: SD of a module-level shift in theta (what the test should detect).
    ``technical``: modules are bins of log_mean_count, and evidence rises with it in BOTH
    halves -- a module that is merely a bin of well-measured genes.
    """
    rng = np.random.default_rng(seed)
    lmc = rng.normal(5, 1, n_genes)
    if technical:
        module = pd.qcut(lmc, n_modules, labels=False)
    else:
        module = rng.integers(0, n_modules, n_genes)
    shift = rng.normal(0, 1, n_modules)[module] * module_signal
    theta = np.clip(0.5 + shift + technical * (lmc - 5) + rng.normal(0, 0.7, n_genes), 0, None)
    e_disc = theta + rng.exponential(0.8, n_genes)
    e_rep = theta + rng.exponential(0.8, n_genes)
    d = pd.DataFrame({
        "gene": [f"G{i}" for i in range(n_genes)], "module_id": module.astype(str),
        "e_disc": e_disc, "e_rep": e_rep, "n_tx": rng.integers(2, 8, n_genes),
        "log_mean_count": lmc, "minor_usage": rng.uniform(0, 0.5, n_genes),
        "p_rep": 10 ** -e_rep, "q_disc": m.bh(10 ** -e_disc),
        "ctx_abund": rng.normal(size=n_genes), "ctx_net": rng.normal(size=n_genes),
    })
    d["ctx_module"] = m.loo_context(d["module_id"].to_numpy(), d["e_disc"].to_numpy())
    d.attrs["n_universe"] = n_genes
    return d


def test_planted_module_signal_is_detected():
    d = _sim(module_signal=0.6, seed=4)
    obs, null = m.replicate_stats(d, n_perm=200, rng=np.random.default_rng(5))
    assert obs["r_module"] > 0.05
    assert m.perm_p(obs["r_module"], null["r_module_strat"]) < 0.01


def test_technical_binning_does_not_pass_the_stratified_null():
    """Modules that are expression bins, with evidence tracking expression, are not signal."""
    d = _sim(module_signal=0.0, technical=0.8, seed=6)
    obs, null = m.replicate_stats(d, n_perm=200, rng=np.random.default_rng(7))
    assert m.perm_p(obs["r_module"], null["r_module_strat"]) > 0.05


def test_no_signal_is_not_detected():
    d = _sim(module_signal=0.0, seed=8)
    obs, null = m.replicate_stats(d, n_perm=200, rng=np.random.default_rng(9))
    assert m.perm_p(obs["r_module"], null["r_module_strat"]) > 0.05


# --------------------------------------------------------------------------- #
# decision rule
# --------------------------------------------------------------------------- #
def _summary(qs):
    cohorts = ["brainseq"] * 3 + ["gtex"] * 3
    return pd.DataFrame({"cohort": cohorts, "q_strat": qs})


def test_verdict_rule():
    assert m.verdict(_summary([.01, .01, .01, .01, .5, .5])) == "SUPPORTED"
    assert m.verdict(_summary([.01, .01, .01, .5, .5, .5])) == "PARTIAL"
    assert m.verdict(_summary([.5] * 5 + [.01])) == "PARTIAL"
    assert m.verdict(_summary([.5] * 6)) == "NOT SUPPORTED"


def test_bh_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = np.array([0.01, 0.04, 0.03, 0.5, 0.2, 0.001])
    assert m.bh(p) == pytest.approx(multipletests(p, method="fdr_bh")[1])


# --------------------------------------------------------------------------- #
# reporting: split distribution, sub-threshold keys, giant-module supplement
# --------------------------------------------------------------------------- #
def test_replicate_stats_reports_subthreshold_enrichment():
    d = _sim(module_signal=0.6, seed=10)
    obs, _ = m.replicate_stats(d, n_perm=10, rng=np.random.default_rng(11))
    for k in ("n_subthreshold", "subthr_rate_top", "subthr_rate_bottom", "subthr_or_top_vs_bottom"):
        assert k in obs
    assert obs["n_subthreshold"] + obs["n_disc_sig"] == len(d)


def test_split_distribution():
    got = m.split_distribution(pd.Series([-0.2, 0.1, 0.1, 0.3]), "r")
    assert got["r"] == pytest.approx(0.075)
    assert got["r_median"] == pytest.approx(0.1)
    assert (got["r_min"], got["r_max"]) == (-0.2, 0.3)
    assert got["r_n_positive"] == 3


def test_region_summary_permutation_statistic_is_the_split_mean():
    rep = pd.DataFrame({
        "cohort": "gtex", "region": "x", "n_genes": 100, "n_universe": 200, "n_modules": 5,
        "r_module": [0.1, 0.3], "r_module_given_abund": [0.1, 0.2], "r_abund": [0.0, 0.1],
        "r_net": 0.0, "r_module_given_net": 0.0, "beta_per_sd": 0.0,
        "seed": [0, 0], "discovery": ["A", "B"]})
    null = pd.DataFrame({"perm": [0, 0, 1, 1], "r_module_strat": [0.0, 0.1, 0.5, 0.5],
                         "r_module_plain": 0.0, "r_module_given_abund_strat": 0.0})
    row = m.region_summary(rep, null)
    assert row["r_module"] == pytest.approx(0.2)
    # per-permutation means are 0.05 and 0.5; only the latter reaches 0.2
    assert row["p_strat"] == pytest.approx(2 / 3)
    assert row["r_module_n_positive"] == 2


def test_region_class():
    assert m.region_class(0.01, 0.01) == "robust to co-expression"
    assert m.region_class(0.01, 0.2) == "not independent of co-expression"
    assert m.region_class(0.2, 0.01) == "not detected"


def test_giant_sensitivity_row_matches_observed_statistic():
    d = _sim(module_signal=0.6, seed=12)
    obs, _ = m.replicate_stats(d, n_perm=5, rng=np.random.default_rng(13))
    part = pd.DataFrame({"gene_id": d["gene"], "module_id": d["module_id"]})
    big = part["module_id"].iloc[0]
    part = pd.concat([part, pd.DataFrame({"gene_id": [f"X{i}" for i in range(1000)],
                                          "module_id": big})], ignore_index=True)
    row = m.giant_sensitivity_row(d, part)
    assert row["r_all"] == pytest.approx(obs["r_module"])
    assert row["n_giant_modules"] == 1 and row["largest_module_genes"] >= 1000
    manual = m.module_context_r(d[d["module_id"] != big])
    assert row["r_no_giant"] == pytest.approx(manual)
    assert row["r_no_largest"] == pytest.approx(manual)
