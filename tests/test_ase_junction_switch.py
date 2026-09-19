"""Tests for the allele-aware junction recount of switch pairs (PI item 10a)."""
from __future__ import annotations

import pandas as pd
import pytest

from isograph_benchmark.real_data import ase_junction_switch as ajs

pysam = pytest.importorskip("pysam")


# --------------------------------------------------------------------------- #
# Junctions are two-sided where exonic sequence is not
# --------------------------------------------------------------------------- #
def test_alt_3prime_site_is_two_sided_by_junction():
    """T1's sequence is a subset of T2's (T2 extends exon 2 upstream), so exonic-unique
    sequence exists on T2's side only -- the one-sided failure of the phASER screen. The
    junctions differ on BOTH sides."""
    i1 = ajs.transcript_introns([(100, 200), (300, 400)])
    i2 = ajs.transcript_introns([(100, 200), (250, 400)])
    s1, s2 = ajs.specific_junctions(i1, i2)
    assert s1 == {(201, 299)} and s2 == {(201, 249)}


def test_exon_skipping_is_two_sided_by_junction():
    i1 = ajs.transcript_introns([(100, 200), (300, 400), (500, 600)])
    i2 = ajs.transcript_introns([(100, 200), (500, 600)])
    s1, s2 = ajs.specific_junctions(i1, i2)
    assert s1 == {(201, 299), (401, 499)} and s2 == {(201, 499)}


def test_pure_truncation_stays_one_sided():
    """Same exons, one transcript shorter: no junction belongs to the short one alone."""
    i1 = ajs.transcript_introns([(100, 200), (300, 400), (600, 700)])
    i2 = ajs.transcript_introns([(100, 200), (300, 400)])
    s1, s2 = ajs.specific_junctions(i1, i2)
    assert s1 == {(401, 599)} and s2 == set()


def test_read_junctions_coordinates_and_anchor():
    # 0-based start 170 -> first block covers 171-200, intron 201-299, then 300-329
    assert ajs.read_junctions([(0, 30), (3, 99), (0, 30)], 170) == [(201, 299)]
    # a 5-base overhang is below the default anchor and is not trusted
    assert ajs.read_junctions([(4, 25), (0, 5), (3, 99), (0, 30)], 190) == []
    assert ajs.read_junctions([(0, 5), (3, 99), (0, 30)], 195, min_anchor=5) == [(201, 299)]


def test_assign_isoform_rules():
    pm = ajs.PairModel("T1|T2", (100, 400),
                       (frozenset({(201, 299)}), frozenset({(201, 249)})),
                       (frozenset({(201, 299)}), frozenset({(201, 249)})))
    assert ajs.assign_isoform(frozenset({(201, 299)}), pm) == 1
    assert ajs.assign_isoform(frozenset({(201, 249)}), pm) == 2
    assert ajs.assign_isoform(frozenset(), pm) is None
    assert ajs.assign_isoform(frozenset({(201, 299), (201, 249)}), pm) == "conflict"
    # a junction the assigned isoform does not splice, inside the pair's span
    assert ajs.assign_isoform(frozenset({(201, 299), (321, 340)}), pm) == "incompatible"
    # ... but one outside the span (another gene) does not disqualify
    assert ajs.assign_isoform(frozenset({(201, 299), (901, 950)}), pm) == 1


def test_merge_mates_and_fragment_hap():
    ok = ajs.ReadRec((), {180: 1}, 1)
    assert ajs.merge_mates([ok, ajs.ReadRec((), {}, 2)])[0] == "wasp_fail"
    assert ajs.merge_mates([ajs.ReadRec((), {180: 1}, None)])[0] == "wasp_missing"
    # mates disagree at a shared site -> the site is dropped, not the fragment
    st, _, sites = ajs.merge_mates([ok, ajs.ReadRec((), {180: 2}, 1)])
    assert st == "ok" and sites == {}
    assert ajs.fragment_hap({}) == 0
    assert ajs.fragment_hap({180: 2, 310: 2}) == 2
    assert ajs.fragment_hap({180: 1, 310: 2}) == "conflict"


# --------------------------------------------------------------------------- #
# End to end on a tiny BAM: every fragment fate once
# --------------------------------------------------------------------------- #
def _read(name, start0, cigar, bases: dict[int, str], vw, *, flag=0, mate=None):
    """bases: 1-based ref position -> base; everything else is 'A'."""
    a = pysam.AlignedSegment()
    a.query_name = name
    a.reference_id = 0
    a.reference_start = start0
    a.cigartuples = cigar
    qlen = sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
    seq = ["A"] * qlen
    q = 0
    r = start0
    for op, n in cigar:
        if op in (0, 7, 8):
            for k in range(n):
                if r + k + 1 in bases:
                    seq[q + k] = bases[r + k + 1]
            q += n
            r += n
        elif op in (2, 3):
            r += n
        elif op in (1, 4):
            q += n
    a.query_sequence = "".join(seq)
    a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    a.mapping_quality = 255
    a.flag = flag
    tags = [("NH", 1)]
    if vw is not None:
        tags.append(("vW", vw))
    a.set_tags(tags)
    if mate is not None:
        a.next_reference_id = 0
        a.next_reference_start = mate
    return a


def test_count_bam_end_to_end(tmp_path):
    t1 = [(0, 30), (3, 99), (0, 30)]      # 171-200 ^ 201-299 ^ 300-329  (T1 junction)
    t2 = [(0, 30), (3, 49), (0, 30)]      # 171-200 ^ 201-249 ^ 250-279  (T2 junction)
    P1, P2 = 0x1 | 0x40, 0x1 | 0x80
    reads = [
        # frag1: T1 junction, hap1 at 180; mate has no site and no WASP tag
        _read("f1", 170, t1, {180: "A"}, 1, flag=P1, mate=350),
        _read("f1", 350, [(0, 50)], {}, None, flag=P2, mate=170),
        # frag2: single read, T2 junction, hap2 at 180
        _read("f2", 170, t2, {180: "G"}, 1),
        # frag3: WASP-fail mate drops the whole fragment
        _read("f3", 170, t1, {180: "A"}, 2),
        # frag4: mates disagree at 180 -> site dropped -> hap 0, still isoform 1
        _read("f4", 170, t1, {180: "A"}, 1, flag=P1, mate=160),
        _read("f4", 160, [(0, 50)], {180: "G"}, 1, flag=P2, mate=170),
        # frag5: an extra junction T1 does not splice -> incompatible
        _read("f5", 170, [(0, 30), (3, 99), (0, 21), (3, 20), (0, 20)], {180: "A"}, 1),
        # frag6: two sites naming different haplotypes -> hap conflict
        _read("f6", 170, t1, {180: "A", 310: "C"}, 1),
        # frag7: covers a het site with no WASP tag -> wasp_missing
        _read("f7", 170, t2, {180: "G"}, None),
    ]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 5000}]}
    unsorted = tmp_path / "u.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for r in sorted(reads, key=lambda x: x.reference_start):
            out.write(r)
    bam = tmp_path / "t.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))

    het = {"chr1": ajs.HetSites(pos=[180, 310],
                                hap_alleles={180: ("A", "G"), 310: ("T", "C")},
                                ref_alt={180: ("A", "G"), 310: ("C", "T")})}
    pm = ajs.PairModel("T1|T2", (100, 400),
                       (frozenset({(201, 299)}), frozenset({(201, 249)})),
                       (frozenset({(201, 299)}), frozenset({(201, 249)})))
    index = {("chr1", 201, 299): [0], ("chr1", 201, 249): [0]}
    pc, ps, sc, qc = ajs.count_bam(str(bam), [("chr1", 0, 1000)], het, [pm], index)

    assert pc == {("T1|T2", 1, 1): 1, ("T1|T2", 2, 2): 1, ("T1|T2", 1, 0): 1}
    assert ps == {("T1|T2", 1, 1, 180): 1, ("T1|T2", 2, 2, 180): 1}
    assert qc["wasp_fail"] == 1
    assert qc["wasp_missing"] == 1
    assert qc["iso_incompatible"] == 1
    assert qc["hap_conflict"] == 1
    assert qc["fragments"] == 7
    # isoform-agnostic site counts: f1 (hap1), f2 (hap2), f5 (hap1); f6 is a conflict
    assert sc == {("chr1", 180, 1): 2, ("chr1", 180, 2): 1}


def test_load_het_sites_uses_pw_phase(tmp_path):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=5000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
        '##FORMAT=<ID=PW,Number=1,Type=String,Description="pw">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBr1\n"
        "chr1\t180\t.\tA\tG\t.\t.\t.\tGT:PW\t0|1:1|0\n"      # PW flips the GT phase
        "chr1\t200\t.\tA\tG\t.\t.\t.\tGT:PW\t0|0:0|0\n"      # homozygous: not a het site
        "chr1\t220\t.\tAT\tG\t.\t.\t.\tGT:PW\t0|1:0|1\n"     # indel: skipped
        "chr1\t240\t.\tC\tT\t.\t.\t.\tGT:PW\t1|0:.\n")       # no PW: GT fallback
    gz = pysam.tabix_index(str(vcf), preset="vcf", force=True)
    sites, qc = ajs.load_het_sites(gz, [("chr1", 0, 5000)])
    h = sites["chr1"]
    assert h.pos == [180, 240]
    assert h.hap_alleles[180] == ("G", "A")
    assert h.hap_alleles[240] == ("T", "C")
    assert qc["het_gt_fallback"] == 1


# --------------------------------------------------------------------------- #
# The gate is the exonic screen's gate
# --------------------------------------------------------------------------- #
def test_gate_requires_both_isoforms_and_donor_depth():
    rows = []
    for d in range(30):
        rows += [("both", f"D{d}", 1, 1, 3), ("both", f"D{d}", 2, 2, 2)]
        rows += [("one_side", f"D{d}", 1, 1, 9)]
        rows += [("shallow", f"D{d}", 1, 1, 2), ("shallow", f"D{d}", 2, 2, 2)]
        rows += [("unphased", f"D{d}", 1, 0, 50), ("unphased", f"D{d}", 2, 0, 50)]
    counts = pd.DataFrame(rows, columns=["pair_id", "donor_id", "isoform", "hap", "n_frag"])
    g = ajs.gate(counts, 30, 5).set_index("pair_id")
    assert bool(g.loc["both", "testable"])
    assert g.loc["both", "n_donors_both_isoforms"] == 30
    assert not bool(g.loc["one_side", "testable"])
    assert not bool(g.loc["shallow", "testable"])
    assert "unphased" not in g.index   # hap 0 carries no allelic information
