"""The BrainSEQ switch-QTL layer rests on one claim that has to be true numerically.

The validation plan asserted 452-500 donors for a BrainSEQ switch-QTL. That figure is
the *junction* phenotype's n; the IsoGraph switch coordinate exists only for the samples
the discovery VAE was fit on -- controls, adults, 238/222/238. The whole expanded-cohort
design depends on `S_g` and `A_g` being **deterministic per-gene transforms of the
transcript counts**, not VAE outputs, so that they can be recomputed on any sample set.

If that were false -- if the switch coordinate carried anything learned from the fit --
then recomputing it on a different sample set would silently change the phenotype and the
QTL would be mapped against something that is not the published switch axis.

These tests pin that property, plus the two places the matched design could quietly
become unmatched (covariates, and the inverse normal transform).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data.brainseq_switch_qtl import (
    AGE_INCLUSIVE,
    AGE_MIN,
    MODALITIES,
    QTL_COVARIATES,
    _covariate_frame,
    _hidden_factors,
    _inverse_normal_transform,
)


def _counts(n_tx: int, n_samples: int, seed: int = 0) -> tuple[np.ndarray, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    m = rng.lognormal(mean=3.0, sigma=1.0, size=(n_tx, n_samples))
    table = pd.DataFrame({
        "transcript_id": [f"ENST{i:08d}" for i in range(n_tx)],
        "gene_id": [f"ENSG{i // 3:08d}" for i in range(n_tx)],
    })
    return m, table


def test_switch_coordinate_is_a_deterministic_function_of_the_counts():
    """Same counts in, bit-identical S_g out -- no hidden state, no RNG."""
    from isograph.features.channels import gene_feature_channels

    m, table = _counts(30, 40)
    a, info_a = gene_feature_channels(m, table)
    b, info_b = gene_feature_channels(m, table)
    np.testing.assert_array_equal(a, b)
    pd.testing.assert_frame_equal(info_a, info_b)


def test_switch_coordinate_of_a_donor_subset_matches_a_refit_on_that_subset():
    """Recomputing on a sample subset is a *refit*, not a projection -- and must be
    reproducible as such.

    This is the property the expanded-cohort design needs: `gene_feature_channels` run on
    a given sample set depends on nothing but that sample set's counts. Running it twice
    on the same subset must agree exactly, which is what makes the QTL phenotype
    regenerable from the committed CLI rather than from a saved model.
    """
    from isograph.features.channels import gene_feature_channels

    m, table = _counts(30, 40, seed=1)
    idx = np.arange(0, 40, 2)
    first, _ = gene_feature_channels(m[:, idx], table)
    second, _ = gene_feature_channels(m[:, idx], table)
    np.testing.assert_allclose(first, second, rtol=0, atol=0)
    # And it is genuinely a refit: the subset's scores are NOT a column slice of the
    # full-cohort scores, because PC1 and the z-scaling are re-estimated on the subset.
    full, _ = gene_feature_channels(m, table)
    assert not np.allclose(first, full[:, idx])


def test_both_modalities_are_produced_for_the_same_genes():
    """A gene with >= 2 transcripts must yield both an S_g and an A_g row.

    The paired swQTL-vs-eQTL test is within-gene, so a gene present on only one axis
    would silently drop out of the contrast rather than count as discordant.
    """
    from isograph.features.channels import gene_feature_channels

    m, table = _counts(30, 40, seed=2)
    _, info = gene_feature_channels(m, table)
    sw = set(info.loc[info["feature_type"] == "switch", "gene_id"])
    ab = set(info.loc[info["feature_type"] == "abundance", "gene_id"])
    assert sw, "no switch features produced"
    assert sw <= ab, "every switch gene must also carry an abundance feature"
    assert set(MODALITIES) == {"switch", "abundance"}


def test_inverse_normal_transform_is_identical_across_modalities():
    """INT is applied to S_g and A_g alike, so a QTL-yield difference between the axes
    cannot be an artefact of one being heavier-tailed than the other."""
    rng = np.random.default_rng(3)
    heavy = pd.DataFrame(rng.standard_cauchy((5, 60)))     # heavy-tailed, like a PC1
    light = pd.DataFrame(rng.normal(size=(5, 60)))
    a = _inverse_normal_transform(heavy)
    b = _inverse_normal_transform(light)
    # After INT both have the same marginal distribution, by construction.
    np.testing.assert_allclose(np.sort(a.to_numpy(), axis=1),
                               np.sort(b.to_numpy(), axis=1), rtol=1e-10)
    assert np.isfinite(a.to_numpy()).all()


def test_inverse_normal_transform_is_rank_preserving():
    rng = np.random.default_rng(4)
    d = pd.DataFrame(rng.normal(size=(3, 50)))
    t = _inverse_normal_transform(d)
    for i in range(len(d)):
        assert (d.iloc[i].rank().to_numpy() == t.iloc[i].rank().to_numpy()).all()


def test_covariate_frame_drops_constant_columns_and_keys_on_donor():
    """A constant covariate makes the design rank-deficient; it must be dropped loudly
    rather than reaching tensorQTL."""
    s = pd.DataFrame({
        "BrNum": [f"Br{i}" for i in range(10)],
        "SNP_PC1": np.linspace(-1, 1, 10), "SNP_PC2": np.linspace(1, -1, 10),
        "SNP_PC3": np.zeros(10),                       # constant -> dropped
        "SNP_PC4": np.linspace(0, 1, 10), "SNP_PC5": np.linspace(1, 2, 10),
        "Sex": ["M", "F"] * 5,
        "Dx": ["Control"] * 5 + ["SCZD"] * 5,
        "RIN": np.linspace(5, 9, 10),
        "mito_rate": np.linspace(0.01, 0.2, 10),
        "mapping_rate": np.linspace(0.8, 0.99, 10),
    })
    cov = _covariate_frame(s)
    assert list(cov.columns) == list(s["BrNum"])
    assert "SNP_PC3" not in cov.index
    assert cov.notna().all().all()
    # categoricals dummied with one level dropped, so the design keeps full rank
    assert "Sex_M" in cov.index and "Sex_F" not in cov.index


def test_covariate_frame_rejects_a_missing_covariate():
    s = pd.DataFrame({"BrNum": ["Br1", "Br2"], "SNP_PC1": [0.0, 1.0]})
    with pytest.raises(SystemExit):
        _covariate_frame(s)


def test_hidden_factors_are_orthogonal_to_the_known_covariates():
    """Factors are PCs of the phenotype matrix AFTER residualizing on the known
    covariates, so they cannot simply re-encode sex, ancestry or RIN and go collinear
    with the covariates already in the model."""
    rng = np.random.default_rng(5)
    donors = [f"Br{i}" for i in range(60)]
    known = pd.DataFrame({"RIN": rng.normal(size=60), "PC1": rng.normal(size=60)},
                         index=donors)
    # phenotypes driven partly by RIN, plus latent structure
    latent = rng.normal(size=60)
    Y = (np.outer(rng.normal(size=200), known["RIN"].to_numpy())
         + np.outer(rng.normal(size=200), latent)
         + rng.normal(size=(200, 60)) * 0.1)
    pheno = pd.DataFrame(Y, columns=donors)
    F = _hidden_factors(pheno, known, k=3)
    assert F.shape == (60, 3)
    for c in F.columns:
        for k in known.columns:
            r = np.corrcoef(F[c], known[k])[0, 1]
            assert abs(r) < 0.15, f"{c} still correlated with {k}: r={r:.3f}"


def test_the_agreed_sample_filter_is_what_the_code_applies():
    """`Age > 13` and `dropped != 't'` were agreed explicitly and are looser than the
    discovery bundle's `Age >= 18, Dx == Control`. Pin them so a later edit to the
    bundle builder cannot silently drag the QTL cohort back to the discovery one."""
    assert AGE_MIN == 13.0
    assert AGE_INCLUSIVE is False
    assert "Dx" in QTL_COVARIATES, "expanded cohort mixes diagnoses; Dx must be adjusted"
    for pc in ("SNP_PC1", "SNP_PC2", "SNP_PC3", "SNP_PC4", "SNP_PC5"):
        assert pc in QTL_COVARIATES, "admixed cohort needs genotype PCs"


def test_a_sub_panel_arm_uses_its_own_genotype_pcs():
    """EA-only mapping adjusts for PCs computed within the EA panel, not the multi-ancestry
    PCs the bundle carries; all_samples keeps the bundle PCs untouched."""
    from isograph_benchmark.real_data.brainseq_switch_qtl import arm_genotype_pcs
    samples = pd.DataFrame({"BrNum": ["Br1", "Br2"], "RIN": [7.0, 8.0],
                            "SNP_PC1": [0.5, -0.5], "SNP_PC2": [0.1, 0.2]})
    ea = pd.DataFrame({"#IID": ["Br2", "Br1", "Br9"], "PC1": [0.02, 0.01, 0.0],
                       "PC2": [0.3, 0.4, 0.0]})
    assert arm_genotype_pcs(samples, "all_samples") is samples
    out = arm_genotype_pcs(samples, "ea_only", pcs=ea).set_index("BrNum")
    assert out.loc["Br1", "SNP_PC1"] == 0.01 and out.loc["Br2", "SNP_PC2"] == 0.3
    assert out.loc["Br1", "RIN"] == 7.0


def test_a_donor_without_sub_panel_pcs_is_refused_not_imputed():
    from isograph_benchmark.real_data.brainseq_switch_qtl import arm_genotype_pcs
    samples = pd.DataFrame({"BrNum": ["Br1", "Br2"], "SNP_PC1": [0.5, -0.5]})
    ea = pd.DataFrame({"#IID": ["Br1"], "PC1": [0.01]})
    with pytest.raises(SystemExit):
        arm_genotype_pcs(samples, "ea_only", pcs=ea)
