"""Tests for the allele-specific switch test on the junction recount (PI item 10a, step 4)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data import ase_junction_allelic as aja

pytest.importorskip("scipy")


# --------------------------------------------------------------------------- #
# Orientation: the lead's phased GT names the ALT haplotype in the PW frame
# --------------------------------------------------------------------------- #
def test_orient_reads_alt_haplotype_off_the_phased_lead_gt():
    rows = []
    for donor in ("D_10", "D_01", "D_00"):
        # T1 on hap 1: 5, T1 on hap 2: 1, T2 on hap 1: 2, T2 on hap 2: 7; hap 0 is ignored
        rows += [("P", donor, 1, 1, 5), ("P", donor, 1, 2, 1),
                 ("P", donor, 2, 1, 2), ("P", donor, 2, 2, 7), ("P", donor, 1, 0, 99)]
    counts = pd.DataFrame(rows, columns=["pair_id", "donor_id", "isoform", "hap", "n_frag"])
    lead = pd.DataFrame({"variant_id": "rs1", "donor_id": ["D_10", "D_01", "D_00"],
                         "gt": ["1|0", "0|1", "0|0"]})
    pairs = pd.DataFrame({"pair_id": ["P"], "gene_id": ["G"], "variant_id_all": ["rs1"]})
    o = aja.orient(counts, lead, pairs).set_index("donor_id")
    # 1|0: ALT is the left allele = hap 1
    assert tuple(o.loc["D_10", ["t1_alt", "t2_alt", "t1_ref", "t2_ref"]]) == (5, 2, 1, 7)
    assert o.loc["D_10", "alt_hap"] == 1
    # 0|1: ALT is hap 2, so the haplotypes swap
    assert tuple(o.loc["D_01", ["t1_alt", "t2_alt", "t1_ref", "t2_ref"]]) == (1, 7, 5, 2)
    assert o.loc["D_01", "alt_hap"] == 2
    # homozygous: no ALT haplotype, hap 1 is the pseudo-ALT of the null
    assert o.loc["D_00", "alt_hap"] == 0
    assert tuple(o.loc["D_00", ["t1_alt", "t1_ref"]]) == (5, 1)


# --------------------------------------------------------------------------- #
# The model recovers a known within-donor effect and is null without one
# --------------------------------------------------------------------------- #
def _simulate(beta, n_donors, mean_n, seed, sigma=1.0, rho=0.05):
    rng = np.random.default_rng(seed)
    rows = []
    for d in range(n_donors):
        u = rng.normal(0, sigma)
        rec = {"donor_id": f"D{d}", "alt_hap": 1}
        for hap, x in (("alt", 1), ("ref", 0)):
            p = 1 / (1 + np.exp(-(0.3 + u + beta * x)))
            n = rng.poisson(mean_n)
            s = (1 - rho) / rho
            q = rng.beta(p * s, (1 - p) * s)
            y = rng.binomial(n, q)
            rec[f"t1_{hap}"], rec[f"t2_{hap}"] = y, n - y
        rows.append(rec)
    return pd.DataFrame(rows)


def test_glmm_recovers_effect():
    tab = _simulate(beta=1.0, n_donors=120, mean_n=10, seed=1)
    f = aja.fit_bb_glmm(*aja._units(tab))
    assert abs(f["beta"] - 1.0) < 0.3
    assert f["pval"] < 1e-6 and f["converged"]


def test_glmm_is_not_biased_away_from_zero_on_sparse_counts():
    """Two units per donor with ~3 fragments each: the regime where fixed donor intercepts
    inflate beta (toward 2x). The random intercept should not."""
    est = [aja.fit_bb_glmm(*aja._units(_simulate(0.8, 300, 3, seed)))["beta"]
           for seed in range(3)]
    assert 0.55 < float(np.mean(est)) < 1.05


def test_glmm_null():
    tab = _simulate(beta=0.0, n_donors=120, mean_n=10, seed=2)
    f = aja.fit_bb_glmm(*aja._units(tab))
    assert abs(f["beta"]) < 0.3 and f["pval"] > 0.01


def test_robust_score_sign_and_null():
    pos = aja.robust_score(_simulate(1.0, 120, 10, seed=3))
    assert pos["z"] > 4
    neg = aja.robust_score(_simulate(-1.0, 120, 10, seed=3))
    assert neg["z"] < -4


# --------------------------------------------------------------------------- #
# Per-pair bookkeeping
# --------------------------------------------------------------------------- #
def test_allelic_contrast_needs_paired_het_donors():
    tab = _simulate(1.0, 8, 10, seed=4)
    out = aja.allelic_contrast(tab, min_paired=10)
    assert out["status"] == "too_few_paired_het_donors"
    assert "beta" not in out
    tab = pd.concat([_simulate(1.0, 40, 10, seed=5),
                     _simulate(0.0, 40, 10, seed=6).assign(alt_hap=0,
                                                          donor_id=lambda d: "H" + d.donor_id)])
    out = aja.allelic_contrast(tab, min_paired=10)
    assert out["status"] == "fitted" and out["n_paired_het"] == 40 and out["n_donors_hom"] == 40
    assert out["beta"] > 0 and "pval_hom_null" in out


def test_between_donor_uses_dosage():
    rng = np.random.default_rng(7)
    dos = pd.Series(rng.integers(0, 3, 60), index=[f"D{i}" for i in range(60)])
    n = 40
    t1 = rng.binomial(n, 0.2 + 0.25 * dos.to_numpy())
    iso = pd.DataFrame({"donor_id": dos.index, "t1_all": t1, "t2_all": n - t1})
    r = aja.between_donor(iso, dos, min_reads=5)
    assert r["rho_between"] > 0.5 and r["pval_between"] < 1e-4


# --------------------------------------------------------------------------- #
# IsoGraph modules
# --------------------------------------------------------------------------- #
def test_add_module_columns_signs_beta_by_module_polarity():
    out = pd.DataFrame({"gene_id": ["G1", "G1", "G2"], "transcript_id_1": ["A", "B", "X"],
                        "transcript_id_2": ["B", "A", "Y"], "beta": [0.5, -0.5, 1.0],
                        "beta_risk": [np.nan, np.nan, -1.0]})
    roles = pd.DataFrame({"gene_id": ["G1"], "module_id": ["M001"], "module_role": ["discordant"]})
    pol = pd.Series({("M001", "A"): 0.6, ("M001", "B"): -0.2})
    age = pd.DataFrame({"module_id": ["M001"], "module_age_trait": ["Age_linear"],
                        "module_age_effect": [0.3]})
    o = aja.add_module_columns(out, roles, pol, age)
    assert o["module_polarity"].round(6).tolist()[:2] == [0.8, -0.8]
    # flipping T1/T2 flips both beta and polarity: the along-module effect is unchanged
    assert o["beta_along_module"].tolist()[:2] == [0.5, 0.5]
    assert pd.isna(o.loc[2, "module_id"])
    assert pd.isna(o.loc[2, "beta_along_module"])


def _direction_frame(aligned: bool, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(40):
        pol = rng.normal(size=5)
        beta = (1.5 * pol + rng.normal(scale=0.3, size=5)) if aligned else rng.normal(size=5)
        for k in range(5):
            rows.append((f"G{g}", beta[k], pol[k], 0.01 if abs(beta[k]) > 1 else 0.5))
    return pd.DataFrame(rows, columns=["gene_id", "beta", "module_polarity", "qval"])


def test_module_direction_check_detects_alignment_and_is_null_otherwise():
    hit = aja.module_direction_check(_direction_frame(True), n_perm=200)
    assert hit["median_abs_rho"] > 0.8 and hit["perm_p"] < 0.01
    assert hit["n_genes_all_sig_pairs_aligned"] >= 0.8 * hit["n_genes_ge2_sig_pairs"]
    miss = aja.module_direction_check(_direction_frame(False, seed=1), n_perm=200)
    assert miss["perm_p"] > 0.05


def test_add_module_columns_guards_the_projection_with_polarity_q():
    """With q-values available, an unanchored pair is not projected (KLC1, stage 06a 5a)."""
    out = pd.DataFrame({"gene_id": ["G1", "G2"], "transcript_id_1": ["A", "X"],
                        "transcript_id_2": ["B", "Y"], "beta": [0.5, 0.5],
                        "beta_risk": [0.5, 0.5]})
    roles = pd.DataFrame({"gene_id": ["G1", "G2"],
                          "module_id": ["M001", "M002"],
                          "module_role": ["coupled", "abundance_only"]})
    pol = pd.DataFrame(
        {"r": [0.6, -0.2, 0.11, -0.06], "qvalue": [1e-6, 1e-4, 0.17, 0.46]},
        index=pd.MultiIndex.from_tuples(
            [("M001", "A"), ("M001", "B"), ("M002", "X"), ("M002", "Y")],
            names=["module_id", "transcript_id"]))
    age = pd.DataFrame({"module_id": ["M001", "M002"],
                        "module_age_trait": ["Age_linear"] * 2,
                        "module_age_effect": [0.3, 0.3]})
    o = aja.add_module_columns(out, roles, pol, age)
    assert o.loc[0, "module_projection_status"] == "projected"
    assert not pd.isna(o.loc[0, "beta_along_module"])
    # G2 joined its module on abundance: no switch axis to project onto
    assert o.loc[1, "module_projection_status"] == "role_not_switching"
    assert pd.isna(o.loc[1, "beta_along_module"])
    assert pd.isna(o.loc[1, "risk_along_module"])
