"""Sign pin and positive-control logic for the BrainSEQ QTL checks."""
import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data.brainseq_qtl_checks import (
    CONTROL_CONCORDANCE_MIN,
    CONTROL_MIN_STRONG_PAIRS,
    CONTROL_PI1_MIN,
    align_effects,
    arm_check_failures,
    control_verdict,
    gene_replication,
    pin_slopes,
    sign_pin_table,
)
from isograph_benchmark.real_data.brainseq_switch_qtl import (
    _hidden_factors,
    _inverse_normal_transform,
)


def _scores(seed: int = 0, n_genes: int = 4, n_samples: int = 40) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.normal(size=(n_genes, n_samples)),
                        index=[f"ENSG{i}" for i in range(n_genes)],
                        columns=[f"R{j}" for j in range(n_samples)])


def test_sign_pin_detects_a_flip_an_unstable_axis_and_too_few_libraries():
    disc = _scores()
    rng = np.random.default_rng(1)
    rec = disc.copy()
    rec.loc["ENSG1"] = -disc.loc["ENSG1"] + rng.normal(scale=0.05, size=disc.shape[1])
    rec.loc["ENSG2"] = rng.normal(size=disc.shape[1])            # a different axis
    rec = rec.drop(columns=disc.columns[:25])                     # 15 shared libraries left
    t = sign_pin_table(disc, rec, min_shared=10).set_index("gene")
    assert t.loc["ENSG0", "sign"] == 1.0 and not t.loc["ENSG0", "flipped"]
    assert t.loc["ENSG1", "sign"] == -1.0 and t.loc["ENSG1", "flipped"]
    assert t.loc["ENSG2", "axis_unstable"]
    t20 = sign_pin_table(disc, rec, min_shared=20).set_index("gene")
    assert not t20["pinnable"].any() and t20["sign"].isna().all()


def test_pinned_slopes_flip_only_flipped_genes_and_never_guess():
    perm = pd.DataFrame({"phenotype_id": ["ENSG0.3", "ENSG1.1", "ENSG9.2"],
                         "slope": [0.4, 0.4, 0.4], "qval": [0.01, 0.01, 0.01]})
    pin = pd.DataFrame({"gene": ["ENSG0", "ENSG1"], "n_shared": [30, 30], "r": [0.99, -0.98],
                        "sign": [1.0, -1.0], "flipped": [False, True],
                        "axis_unstable": [False, False]})
    p = pin_slopes(perm, pin).set_index("phenotype_id")
    assert p.loc["ENSG0.3", "slope_pinned"] == 0.4
    assert p.loc["ENSG1.1", "slope_pinned"] == -0.4
    assert np.isnan(p.loc["ENSG9.2", "slope_pinned"]) and p.loc["ENSG9.2", "slope"] == 0.4


def test_post_hoc_pinning_is_exact_int_is_odd_and_hidden_factors_ignore_sign():
    rng = np.random.default_rng(3)
    donors = [f"Br{i}" for i in range(50)]
    x = pd.DataFrame(np.round(rng.normal(size=(30, 50)), 1), columns=donors)  # ties on purpose
    assert np.allclose(_inverse_normal_transform(-x), -_inverse_normal_transform(x))
    known = pd.DataFrame({"RIN": rng.normal(size=50)}, index=donors)
    flip = np.where(rng.random(30) < 0.5, -1.0, 1.0)
    F1 = _hidden_factors(x, known, k=3)
    F2 = _hidden_factors(x.mul(flip, axis=0), known, k=3)
    assert np.allclose(np.abs(F1.to_numpy()), np.abs(F2.to_numpy()), atol=1e-8)


def _pairs(**kw):
    base = dict(ref="A", alt="G", bs_ref="A", bs_alt="G", slope_bs=0.5, slope_gtex=0.3)
    base.update(kw)
    return base


def test_alignment_keeps_same_flips_swapped_and_drops_what_it_cannot_resolve():
    pairs = pd.DataFrame([
        _pairs(),                                                   # same orientation
        _pairs(bs_ref="G", bs_alt="A", slope_bs=-0.5),              # swapped -> flip to +
        _pairs(ref="A", alt="T", bs_ref="A", bs_alt="T"),           # palindromic, same: keep
        _pairs(ref="A", alt="T", bs_ref="T", bs_alt="A"),           # palindromic swap: drop
        _pairs(bs_alt="C"),                                         # allele mismatch: drop
        _pairs(slope_bs=-0.5),                                      # discordant
    ])
    out = align_effects(pairs)
    assert list(out.index) == [0, 1, 2, 5]
    assert out.loc[1, "slope_bs_aligned"] == 0.5
    assert out["concordant"].tolist() == [True, True, True, False]


def test_gene_replication_reads_brainseq_signal_at_gtex_egenes():
    rng = np.random.default_rng(4)
    genes = [f"ENSG{i}" for i in range(400)]
    gtex = pd.DataFrame({"gene": genes, "qval": np.r_[np.full(200, 0.001), np.full(200, 0.5)],
                         "pval_nominal": rng.random(400), "pval_beta": rng.random(400)})
    bs = pd.DataFrame({"phenotype_id": [f"{g}.1" for g in genes],
                       "pval_beta": np.r_[rng.random(200) * 1e-4, rng.random(200)],
                       "qval": np.r_[np.full(200, 0.01), np.full(200, 0.9)]})
    r = gene_replication(gtex, bs, top_n=50)
    assert r["n_egenes_tested_in_brainseq"] == 200
    assert r["pi1_brainseq_given_gtex_egene"] > 0.9
    assert r["frac_brainseq_q05_top_gtex_egenes"] == 1.0


def test_control_verdict_applies_the_prespecified_rule():
    assert control_verdict(0.8, 0.97, 500) == ("PASS", "")
    assert control_verdict(CONTROL_PI1_MIN - 0.01, 0.97, 500)[0] == "FAIL"
    assert control_verdict(0.8, CONTROL_CONCORDANCE_MIN - 0.01, 500)[0] == "FAIL"
    v, why = control_verdict(0.8, 1.0, CONTROL_MIN_STRONG_PAIRS - 1)
    assert v == "FAIL" and "pairs" in why
    assert control_verdict(np.nan, 0.97, 500)[0] == "FAIL"


def test_downstream_gate_needs_both_checks_run_and_passed_in_every_region(tmp_path):
    regions = ("caudate", "dlpfc")
    assert len(arm_check_failures("ea_only", regions, tmp_path)) == 2   # neither summary exists

    pd.DataFrame({"region": regions, "arm": "ea_only", "reproduces_discovery": [True, True],
                  "median_abs_r": [0.9996, 0.9996]}).to_parquet(tmp_path / "sign_pin_summary.parquet")
    pd.DataFrame({"region": regions, "arm": "ea_only", "verdict": ["PASS", "FAIL"],
                  "fail_reasons": ["", "pi1 0.30 < 0.5"]}).to_parquet(
        tmp_path / "positive_control_summary.parquet")
    f = arm_check_failures("ea_only", regions, tmp_path)
    assert len(f) == 1 and f[0].startswith("dlpfc") and "pi1" in f[0]

    # A region the checks never covered is a failure, not a pass by omission.
    assert len(arm_check_failures("ea_only", (*regions, "hippocampus"), tmp_path)) == 3
    # Another arm's summary rows do not certify this arm.
    assert len(arm_check_failures("all_samples", regions, tmp_path)) == 4
