"""Structural half of the switch-pair assayability call.

The expression half needs the GTEx TPM tables and is exercised by the SLURM wrapper; what is
worth pinning here is the geometry, because `_classify` decides whether a pair reaches the
bench at all and its failure mode is silent -- a mis-classified pair just disappears from the
panel with a plausible-looking label.
"""
import pytest

from isograph.explain.structure import TranscriptRecord
from isograph_benchmark.real_data.rbp_pair_assayability import _classify, _introns


def _tx(name, exons):
    return TranscriptRecord(transcript_id=name, gene_id="G", chrom="chr1", strand="+",
                            biotype="protein_coding", exons=list(exons), cds=[])


def test_introns_are_the_gaps_between_sorted_exons():
    assert _introns([(300, 400), (100, 200), (500, 600)]) == {(200, 300), (400, 500)}
    assert _introns([(100, 200)]) == set()


def test_cassette_exon_gives_a_two_sided_junction_assay():
    """Exon skipping: the inclusion isoform owns two junctions, the skipping isoform one."""
    inc = _tx("inc", [(100, 200), (300, 400), (500, 600)])
    skip = _tx("skip", [(100, 200), (500, 600)])
    r = _classify(inc, skip, min_unique_bp=80)
    assert r["assay_class"] == "junction"
    assert r["structurally_measurable"]
    assert r["n_unique_junctions_1"] == 2 and r["n_unique_junctions_2"] == 1


def test_retained_intron_is_measurable_despite_one_sided_junctions():
    """The regression this test exists for.

    A retained-intron pair has no junction unique to the retaining isoform -- it has no
    junction at all. Requiring unique exonic bp on *both* sides threw the pair out, even
    though the retained intron body identifies one isoform and the spliced junction the
    other. Each side needs *a* feature, not the same kind of feature.
    """
    retained = _tx("ri", [(100, 600)])
    spliced = _tx("sp", [(100, 200), (500, 600)])
    r = _classify(retained, spliced, min_unique_bp=80)
    assert r["assay_class"] == "junction_and_segment"
    assert r["structurally_measurable"]
    assert r["unique_exonic_bp_1"] == 300      # the intron body, unique to the retaining form
    assert r["unique_exonic_bp_2"] == 0        # spliced form is a subset of the retained one
    assert r["n_unique_junctions_2"] == 1


def test_alternate_terminal_exons_with_unique_sequence_are_segment_assays():
    """Same junction, different ends, both ends long enough to amplify inside."""
    a = _tx("a", [(100, 300), (500, 600)])
    b = _tx("b", [(200, 300), (500, 700)])
    r = _classify(a, b, min_unique_bp=80)
    assert r["assay_class"] == "unique_segment"
    assert r["structurally_measurable"]
    assert r["n_unique_junctions_1"] == 0 and r["n_unique_junctions_2"] == 0


def test_nested_terminal_extension_is_not_measurable_by_a_ratio_assay():
    """A 3' extension leaves the shorter isoform with nothing of its own to amplify."""
    short = _tx("short", [(100, 200), (300, 400)])
    long_ = _tx("long", [(100, 200), (300, 500)])
    r = _classify(short, long_, min_unique_bp=80)
    assert r["assay_class"] == "nested_terminal"
    assert not r["structurally_measurable"]
    assert r["unique_exonic_bp_1"] == 0


def test_unique_stretch_below_the_primer_floor_fails():
    """A 12 bp difference is real but cannot hold a primer, so it is not an assay."""
    a = _tx("a", [(100, 300), (500, 600)])
    b = _tx("b", [(288, 300), (500, 600)])
    r = _classify(a, b, min_unique_bp=80)
    assert r["unique_exonic_bp_1"] == 188 and r["unique_exonic_bp_2"] == 0
    assert not r["structurally_measurable"]


def test_classification_is_symmetric_in_the_two_transcripts():
    """Which transcript the panel happens to list first must not change the verdict."""
    a = _tx("a", [(100, 200), (300, 400), (500, 600)])
    b = _tx("b", [(100, 200), (500, 600)])
    fwd, rev = _classify(a, b, 80), _classify(b, a, 80)
    assert fwd["assay_class"] == rev["assay_class"]
    assert fwd["structurally_measurable"] == rev["structurally_measurable"]
    assert fwd["unique_exonic_bp_1"] == rev["unique_exonic_bp_2"]


@pytest.mark.parametrize("floor", [40, 80, 150])
def test_min_unique_bp_only_ever_removes_pairs(floor):
    """Raising the primer floor must not promote a pair into a measurable class."""
    a = _tx("a", [(100, 300), (500, 600)])
    b = _tx("b", [(200, 300), (500, 700)])
    r = _classify(a, b, min_unique_bp=floor)
    assert r["structurally_measurable"] == (floor <= 100)
