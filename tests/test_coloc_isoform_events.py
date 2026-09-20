"""The signal-level event layer names the same introns the locus event audit tiers on."""
import pandas as pd

from isograph_benchmark.paths import stage_out
from isograph_benchmark.real_data.coloc_isoform_events import (
    _parse_intron_event,
    signal_event_cells,
    signal_events_path,
)

GENE = "ENSG00000145335"
PID_A = f"chr4:89835692:89836127:clu_43552_-:{GENE}.17"
PID_B = f"chr4:89835692:89836743:clu_43552_-:{GENE}.17"
PID_REP = f"chr4:89726660:89729194:clu_43551_-:{GENE}.17"
KEY = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4", gene=GENE)


def _nom() -> pd.DataFrame:
    return pd.DataFrame([{**KEY, "PP4_sQTL": 0.97}])


def _cells() -> pd.DataFrame:
    return pd.DataFrame([
        {**KEY, "symbol": "SNCA", "tissue": "Brain_Cortex", "PP4_sQTL": 0.97,
         "estimator_sQTL": "susie", "phenotype_id": PID_A},
        {**KEY, "symbol": "SNCA", "tissue": "Brain_Putamen_basal_ganglia", "PP4_sQTL": 0.85,
         "estimator_sQTL": "abf", "phenotype_id": None},
        {**KEY, "symbol": "SNCA", "tissue": "Brain_Hippocampus", "PP4_sQTL": 0.40,
         "estimator_sQTL": "susie", "phenotype_id": PID_B},
        {**KEY, "LOCUS_ID": "locus99_chr4", "symbol": "SNCA", "tissue": "Brain_Cortex",
         "PP4_sQTL": 0.99, "estimator_sQTL": "susie", "phenotype_id": PID_B},
    ])


SREP = pd.DataFrame([
    {"gene": GENE, "tissue": "Brain_Putamen_basal_ganglia", "phenotype_id": PID_REP},
    {"gene": GENE, "tissue": "Brain_Cortex", "phenotype_id": PID_REP},
])


def test_selects_called_tissues_at_the_nominated_locus_only():
    ec = signal_event_cells(_nom(), _cells(), SREP)
    assert sorted(ec["tissue"]) == ["Brain_Cortex", "Brain_Putamen_basal_ganglia"]
    assert set(ec["LOCUS_ID"]) == {"locus08_chr4"}


def test_fitted_intron_is_kept_over_the_representative():
    ec = signal_event_cells(_nom(), _cells(), SREP).set_index("tissue")
    assert ec.loc["Brain_Cortex", "phenotype_id"] == PID_A
    assert ec.loc["Brain_Cortex", "intron_event"] == "chr4:89835692-89836127(-)"


def test_abf_fallback_names_the_representative_intron_and_keeps_its_estimator():
    r = signal_event_cells(_nom(), _cells(), SREP).set_index("tissue").loc[
        "Brain_Putamen_basal_ganglia"]
    assert r["phenotype_id"] == PID_REP
    assert r["estimator_sQTL"] == "abf"


def test_intron_event_parses_like_a_clpp_junction():
    for ev in signal_event_cells(_nom(), _cells(), SREP)["intron_event"]:
        chrom, s, e = _parse_intron_event(ev)
        assert chrom == "chr4" and s < e


def test_a_cell_with_no_intron_at_all_is_dropped_not_invented():
    ec = signal_event_cells(_nom(), _cells(), SREP[SREP["tissue"] == "Brain_Cortex"])
    assert "Brain_Putamen_basal_ganglia" not in set(ec["tissue"])


def test_signal_table_never_lands_on_the_clpp_table():
    p = signal_events_path()
    assert p != stage_out("anchoring.coloc", "coloc_isoform_events_combined.parquet")
    assert "all_introns" in p.parts


def test_cross_tissue_exception_takes_pairs_only_from_regions_carrying_the_junction():
    from isograph_benchmark.real_data.coloc_isoform_events import (
        CROSS_TISSUE_EXCEPTIONS,
        cross_tissue_switch_pairs,
    )
    assert "UNC13A" in CROSS_TISSUE_EXCEPTIONS
    unc13a = "ENSG00000130477"
    by_region = {
        "cerebellum": {},
        "frontal_cortex_ba9": {unc13a: {"T1.1", "T2.1", "T3.1"}},
        "hippocampus": {unc13a: {"T1.1", "T4.1"}},
        "caudate_basal_ganglia": {unc13a: {"T5.1", "T6.1"}},
    }
    regions, members = cross_tissue_switch_pairs({"T1.1", "T9.1"}, unc13a, by_region)
    assert regions == ["frontal_cortex_ba9", "hippocampus"]
    assert members == {"T1.1", "T2.1", "T3.1", "T4.1"}
    assert cross_tissue_switch_pairs({"T9.1"}, unc13a, by_region) == ([], set())


def test_structural_flags_survive_the_nullable_boolean_dtype():
    """`structure_annotations` stores flags as pandas `boolean`, not Python bool.

    `itertuples` yields `numpy.bool_` for that dtype, and `numpy.bool_(True) is True`
    is False, so an identity test drops every flag and every event comes back
    "no annotated structural change". Missing values must still be rejected, which
    rules out a bare `bool(v)` because `bool(pd.NA)` raises.
    """
    from isograph_benchmark.real_data.coloc_isoform_events import (
        _flag_true,
        _switch_struct_label,
    )

    a = pd.DataFrame({
        "transcript_id": ["T1.1", "T2.1"],
        "first_exon_changed": pd.array([True, False], dtype="boolean"),
        "last_exon_changed": pd.array([False, False], dtype="boolean"),
        "internal_exon_difference": pd.array([False, True], dtype="boolean"),
        "cds_changed": pd.array([True, pd.NA], dtype="boolean"),
        "utr_changed": pd.array([pd.NA, pd.NA], dtype="boolean"),
        "biotype_switch": pd.array([False, False], dtype="boolean"),
        "coding_status_change": pd.array([False, False], dtype="boolean"),
    })
    row = next(a.itertuples())
    assert row.first_exon_changed is not True       # the shape of the original bug
    assert _flag_true(row.first_exon_changed)
    assert not _flag_true(row.last_exon_changed)
    assert not _flag_true(pd.NA) and not _flag_true(None)

    class _Ev:
        struct = {r.transcript_id: r for r in a.itertuples()}

    label = _switch_struct_label(_Ev(), {"T1.1", "T2.1"})
    assert label == "cds, first exon, internal exon"
    assert _switch_struct_label(_Ev(), set()) == "no annotated structural change"
