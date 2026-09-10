"""The tier rule decides what may be called a recovered mechanism, so it is pinned here.

This module is the difference between "UNC13A colocalizes with a brain sQTL" and "we
recovered the TDP-43 cryptic-exon mechanism". Those are different claims and only one of
them is publishable as a positive control, so the rule that separates them cannot be
allowed to drift silently.

Two gates guard the top tier and they are INDEPENDENT:
  * `status`         -- were the coordinates verified in this repo (provenance)
  * `evidence_class` -- how strong is the underlying literature claim

The second exists because a TWAS association does not establish a shared causal variant:
LD-driven co-regulation at a locus produces the same signal. Colocalization is the more
conservative test, so a TWAS-only anchor may never certify a mechanism however well its
coordinates line up. That asymmetry is the single most load-bearing line in the module and
most of the tests below exist to hold it in place.
"""
from __future__ import annotations

import pandas as pd
import pytest

from isograph_benchmark.real_data.locus_event_audit import (
    PROMOTING_EVIDENCE,
    PROMOTING_STATUS,
    TIERS,
    _parse_gtex_intron,
    assign_tier,
    interval_matches,
    load_curated_events,
)

# The real UNC13A numbers, so the fixtures cannot drift away from the curated registry.
UNC13A_CE = (17641557, 17642844)          # curated, GENCODE-derived
UNC13A_GTEX_CE = (17641556, 17642845)     # the GTEx LeafCutter phenotype at that intron
UNC13A_COLOC = (17630750, 17632782)       # what actually colocalized


def _curated(**kw) -> pd.DataFrame:
    base = dict(gene="UNC13A", trait="als", status="reviewed",
                evidence_class="functional_validation", chrom="chr19",
                start=UNC13A_CE[0], end=UNC13A_CE[1],
                match_mode="exact", match_tolerance=2)
    return pd.DataFrame([{**base, **kw}])


def _introns(*iv, chrom="chr19"):
    return [(chrom, s, e, f"{chrom}:{s}:{e}:clu_x_-") for s, e in iv]


# --------------------------------------------------------------------------- #
# Interval matching
# --------------------------------------------------------------------------- #
def test_the_1bp_leafcutter_offset_is_absorbed_not_rejected():
    """GTEx's phenotype at the cryptic-exon intron is off by exactly 1 bp per end from
    GENCODE. That is a convention difference, not a different event, and the audit must
    match it -- otherwise the all-introns arm could never confirm UNC13A."""
    assert interval_matches(*UNC13A_GTEX_CE, *UNC13A_CE, mode="exact", tol=2)


def test_the_actually_colocalizing_intron_does_not_match():
    """The intron that carries UNC13A's PP4=0.970 is ~9 kb from the cryptic exon."""
    assert not interval_matches(*UNC13A_COLOC, *UNC13A_CE, mode="exact", tol=2)


def test_exact_mode_rejects_an_intron_sharing_only_one_endpoint():
    """Introns in one LeafCutter cluster are ALTERNATIVE splice choices and routinely
    share an endpoint. An overlap rule would let one certify as the other; exact must
    not. This is the real PICALM shape: 85,974,812-85,981,129 vs 85,978,093-85,981,129."""
    assert not interval_matches(85974812, 85981129, 85978093, 85981129,
                                mode="exact", tol=2)
    # ...and the looser mode would have accepted it, which is why exact is the default.
    assert interval_matches(85974812, 85981129, 85978093, 85981129,
                            mode="overlap", tol=2)


def test_exact_is_the_default_mode():
    cur = _curated()
    cur = cur.drop(columns=["match_mode"])
    v = assign_tier(_introns((85974812, 85981129)), cur)
    assert v["curated_interval_matched"] is False


# --------------------------------------------------------------------------- #
# The evidence gate: TWAS may never promote
# --------------------------------------------------------------------------- #
def test_a_twas_anchor_never_promotes_even_on_a_perfect_coordinate_match():
    """The load-bearing rule. TWAS does not establish a shared causal variant."""
    cur = _curated(gene="PICALM", trait="ad", evidence_class="twas_association",
                   chrom="chr11", start=86026368, end=86031468)
    v = assign_tier(_introns((86026367, 86031469), chrom="chr11"), cur)
    assert v["curated_interval_matched"] is True, "the coordinates DO match"
    assert v["curated_matched"] is False, "but it must not certify a mechanism"
    assert v["tier"] == "disease_locus_splice_linked"


def test_functional_and_coloc_evidence_do_promote():
    for ec in ("functional_validation", "coloc_association"):
        v = assign_tier(_introns(UNC13A_GTEX_CE), _curated(evidence_class=ec))
        assert v["curated_matched"] is True, ec
        assert v["tier"] == "known_mechanism_recovered", ec


def test_twas_is_deliberately_absent_from_the_promoting_set():
    assert "twas_association" not in PROMOTING_EVIDENCE
    assert set(PROMOTING_EVIDENCE) == {"functional_validation", "coloc_association"}


# --------------------------------------------------------------------------- #
# The status gate
# --------------------------------------------------------------------------- #
def test_a_proposed_record_never_promotes():
    """`proposed` means the coordinates were not verified here. Even a perfect match and
    the strongest evidence class must not certify on an unverified interval."""
    v = assign_tier(_introns(UNC13A_GTEX_CE),
                    _curated(status="proposed", evidence_class="functional_validation"))
    assert v["curated_interval_matched"] is True
    assert v["curated_matched"] is False
    assert v["tier"] == "disease_locus_splice_linked"
    assert PROMOTING_STATUS == ("reviewed",)


def test_an_unpinned_record_cannot_be_matched_against():
    """A record with NaN coordinates must be skipped, not crash and not match."""
    cur = _curated(status="proposed", start=float("nan"), end=float("nan"))
    v = assign_tier(_introns(UNC13A_GTEX_CE), cur)
    assert v["curated_interval_matched"] is False
    assert v["tier"] == "disease_locus_splice_linked"


# --------------------------------------------------------------------------- #
# Tiering
# --------------------------------------------------------------------------- #
def test_no_curated_event_is_a_novel_candidate_not_a_failure():
    v = assign_tier(_introns((1000, 2000)), pd.DataFrame(columns=_curated().columns))
    assert v["tier"] == "novel_splice_led_candidate"


def test_no_named_intron_is_not_resolved():
    v = assign_tier([], _curated())
    assert v["tier"] == "not_resolved"
    assert v["curated_matched"] is False


def test_a_curated_event_on_another_chromosome_does_not_match():
    v = assign_tier(_introns(UNC13A_GTEX_CE, chrom="chr7"), _curated())
    assert v["curated_interval_matched"] is False


def test_every_tier_returned_is_declared():
    cases = [([], _curated()),
             (_introns(UNC13A_GTEX_CE), _curated()),
             (_introns(UNC13A_COLOC), _curated()),
             (_introns((1, 2)), pd.DataFrame(columns=_curated().columns))]
    for introns, cur in cases:
        assert assign_tier(introns, cur)["tier"] in TIERS


# --------------------------------------------------------------------------- #
# The registry itself
# --------------------------------------------------------------------------- #
def test_the_committed_registry_loads_and_is_internally_consistent():
    ev = load_curated_events()
    assert len(ev) >= 2
    assert {"UNC13A", "PICALM"} <= set(ev["gene"])
    # every reviewed record must carry coordinates -- enforced by the loader
    rev = ev[ev["status"] == "reviewed"]
    assert not rev[["chrom", "start", "end"]].isna().any().any()
    # and every record must declare how strong its evidence is
    assert ev["evidence_class"].notna().all()


def test_picalm_is_recorded_as_twas_and_so_cannot_promote():
    """Guards the specific finding: PICALM's anchor is a TWAS association, so however
    well the audit matches its interval it stays below the top tier."""
    ev = load_curated_events()
    p = ev[ev["gene"] == "PICALM"].iloc[0]
    assert p["evidence_class"] == "twas_association"
    assert p["evidence_class"] not in PROMOTING_EVIDENCE


def test_gtex_intron_id_parsing():
    assert _parse_gtex_intron(
        "chr19:17641556:17642845:clu_28407_-:ENSG00000130477.16"
    ) == ("chr19", 17641556, 17642845)
    assert _parse_gtex_intron("nonsense") is None
    assert _parse_gtex_intron(None) is None
