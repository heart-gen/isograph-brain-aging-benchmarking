"""The orientation is a sign flip, so the tests are about the sign being right."""
import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data import ase_risk_orientation as aro

# plink2 --ld output for two variants in positive phase between the MINOR alleles.
LD_LOG = """
--ld rsLEAD rsGWAS:

rsLEAD alleles:
  MAJOR = REF = C
  MINOR = T
rsGWAS alleles:
  MAJOR = T
  MINOR = C
  (REF = C)

2438 valid samples; 0.984 het pairs statistically phased.

  r^2 = 0.0699037    |D'| = 0.786856

        Frequencies      :             rsGWAS
  (expectations under LE)          MAJOR       MINOR
                                 ----------  ----------
                           MAJOR  0.587240    0.341742
                                 (0.553919)  (0.375062)
                rsLEAD
                           MINOR  0.009026    0.061993
                                 (0.042346)  (0.028673)

  Major alleles are in phase with each other.
"""


def test_parse_ld_keeps_each_variants_own_alleles():
    ld = aro._parse_ld(LD_LOG)
    assert ld["a"] == {"major": "C", "minor": "T"}
    assert ld["b"] == {"major": "T", "minor": "C"}


def test_signed_r_matches_plinks_r2_and_its_phase_statement():
    ld = aro._parse_ld(LD_LOG)
    r = aro.signed_r(ld, "C", "T")           # the two major alleles
    assert r == pytest.approx(0.2644, abs=1e-3)
    assert r ** 2 == pytest.approx(0.0699037, abs=1e-4)   # plink's own r^2
    assert r > 0                                          # "major alleles are in phase"
    # swapping one side's allele flips the sign and nothing else (the printed frequencies
    # are rounded to 6 dp, so the magnitudes agree to ~1e-5, not exactly)
    assert aro.signed_r(ld, "C", "C") == pytest.approx(-r, abs=1e-4)
    assert aro.signed_r(ld, "T", "C") == pytest.approx(r, abs=1e-4)


def test_signed_r_rejects_an_allele_neither_variant_carries():
    assert aro.signed_r(aro._parse_ld(LD_LOG), "G", "T") is None


def _row(**kw):
    base = dict(beta=1.0, module_polarity=1.0, risk_allele="A", r_lead_risk=0.9,
                status="fitted")
    base.update(kw)
    return base


def test_orientation_flips_only_when_the_risk_allele_rides_the_ref_allele():
    t = aro.orient(pd.DataFrame([
        _row(r_lead_risk=0.9),    # risk travels with lead ALT -> beta kept
        _row(r_lead_risk=-0.9),   # risk travels with lead REF -> beta flipped
    ]))
    assert list(t["beta_risk"]) == [1.0, -1.0]
    assert list(t["orientation_status"]) == ["oriented", "oriented"]


def test_weak_ld_leaves_the_pair_unoriented_rather_than_guessing():
    t = aro.orient(pd.DataFrame([_row(r_lead_risk=0.3)]))
    assert np.isnan(t["beta_risk"].iloc[0])
    assert np.isnan(t["risk_along_module"].iloc[0])
    assert t["orientation_status"].iloc[0] == "r_gate"


def test_missing_risk_allele_and_missing_ld_are_reported_apart():
    t = aro.orient(pd.DataFrame([_row(risk_allele=None), _row(r_lead_risk=np.nan)]))
    assert list(t["orientation_status"]) == ["no_gwas_risk_allele", "no_ld"]


def test_module_polarity_sets_the_disease_direction():
    """risk_along_module > 0 means the risk allele pushes the isoforms the module's way."""
    t = aro.orient(pd.DataFrame([
        _row(module_polarity=1.0),    # T1 rises with the module score
        _row(module_polarity=-1.0),   # T1 falls with it: same beta, opposite meaning
    ]))
    assert list(t["risk_along_module"]) == [1.0, -1.0]


def test_palindromic_gwas_lead_is_left_unoriented():
    """A/T and C/G leads match their alleles on either strand, so a strand disagreement
    between the sumstats and the panel would flip the sign without any error."""
    t = aro.orient(pd.DataFrame([
        _row(risk_allele="A", gwas_other_allele="T"),   # palindromic
        _row(risk_allele="A", gwas_other_allele="G"),   # not
    ]))
    assert list(t["orientation_status"]) == ["palindromic_gwas_lead", "oriented"]
    assert np.isnan(t["beta_risk"].iloc[0])
