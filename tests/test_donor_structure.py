"""Donor structure: region is not analysis, donor is not sample, ids are cohort-local.

BrainSEQ caudate is one region carrying two analyses (aging over the controls, SCZD over
those same control samples plus the patients), so the region and analysis levels must not
be conflated -- that is the correction this module exists to encode.
"""
import pandas as pd
import pytest

ds = pytest.importorskip("isograph_benchmark.real_data.donor_structure")

# Two analyses over one region, plus a second region and a colliding id in another cohort.
COUNTS = pd.DataFrame({
    "analysis": ["caudate (aging)", "caudate (SCZD)", "DLPFC", "GTEx cortex"],
    "region": ["BS caudate", "BS caudate", "BS DLPFC", "GTEx cortex"],
    "cohort": ["BrainSEQ", "BrainSEQ", "BrainSEQ", "GTEx"],
    "phenotype": ["aging", "diagnosis", "aging", "aging"],
    "n_samples": [2, 3, 2, 2],
    "n_donors": [2, 3, 2, 2],
})
DONORS = {
    "caudate (aging)": {"Br1", "Br2"},
    "caudate (SCZD)": {"Br1", "Br2", "Br3"},   # controls + one patient
    "DLPFC": {"Br2", "Br4"},
    "GTEx cortex": {"Br1", "G1"},              # "Br1" here is a DIFFERENT person
}
SAMPLES = {
    "caudate (aging)": {"R1", "R2"},
    "caudate (SCZD)": {"R1", "R2", "R3"},      # same RNums re-used, same region
    "DLPFC": {"R4", "R5"},
    "GTEx cortex": {"S1", "S2"},
}
COHORT_R = dict(zip(COUNTS["region"], COUNTS["cohort"]))
COHORT_A = dict(zip(COUNTS["analysis"], COUNTS["cohort"]))


def _region_donors():
    return ds._union_by(COUNTS["region"], COUNTS["analysis"], DONORS)


def test_region_donor_set_is_the_union_over_its_analyses():
    rd = _region_donors()
    assert rd["BS caudate"] == {"Br1", "Br2", "Br3"}
    assert len(rd) == 3  # three regions, not four analyses


def test_sample_ids_may_repeat_within_a_region_but_never_across_regions():
    notes = ds._check_identifier_levels(SAMPLES, COUNTS)
    assert any("share 2 tissue samples" in n and "BS caudate" in n for n in notes)

    bad = dict(SAMPLES, **{"DLPFC": {"R1", "R5"}})  # an RNum leaking into another region
    with pytest.raises(SystemExit, match="shared across regions"):
        ds._check_identifier_levels(bad, COUNTS)


def test_overlap_is_symmetric_and_never_crosses_cohorts():
    ov = ds._overlap(_region_donors(), COHORT_R, "region").set_index(
        ["region_a", "region_b"])

    assert ov.loc[("BS caudate", "BS DLPFC"), "n_shared_donors"] == 1     # Br2
    assert ov.loc[("BS DLPFC", "BS caudate"), "n_shared_donors"] == 1
    assert ov.loc[("BS caudate", "BS DLPFC"), "jaccard"] == pytest.approx(1 / 4)
    # A colliding identifier string across cohorts is two different people.
    assert ov.loc[("BS caudate", "GTEx cortex"), "n_shared_donors"] == 0
    # The diagonal is the region's own donor count.
    assert ov.loc[("BS caudate", "BS caudate"), "n_shared_donors"] == 3


def test_per_donor_separates_regions_from_analyses():
    inc = ds._incidence(_region_donors(), COHORT_R)
    per = ds._per_donor(inc, DONORS, COUNTS).set_index(["cohort", "donor_id"])

    # Br1 gives one region (caudate) but is carried by two analyses.
    assert per.loc[("BrainSEQ", "Br1"), "n_regions"] == 1
    assert per.loc[("BrainSEQ", "Br1"), "n_analyses"] == 2
    # Br2 gives two regions and three analyses.
    assert per.loc[("BrainSEQ", "Br2"), "n_regions"] == 2
    assert per.loc[("BrainSEQ", "Br2"), "n_analyses"] == 3
    # The same string in the other cohort is a separate person.
    assert per.loc[("GTEx", "Br1"), "n_regions"] == 1


def test_incidence_has_one_row_per_donor_and_region():
    inc = ds._incidence(_region_donors(), COHORT_R)
    assert set(inc.columns) == {"cohort", "donor_id", "region"}
    assert not inc.duplicated(subset=["cohort", "donor_id", "region"]).any()
    assert len(inc) == 3 + 2 + 2  # caudate 3, DLPFC 2, GTEx cortex 2
