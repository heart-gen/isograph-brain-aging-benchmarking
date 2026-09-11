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

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data.locus_event_audit import (
    CONTEXT_DISTINCT_ARM,
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
# --------------------------------------------------------------------------- #
# context_distinct_splice_colocalization
# --------------------------------------------------------------------------- #
def test_the_all_introns_arm_reads_the_all_introns_results(monkeypatch):
    """The two arms write identically-named files into different directories. If the arm
    governed the tier while the numbers came from the representative directory, the audit
    would certify event specificity from a run that never tested the curated event -- the
    exact failure the tier exists to rule out."""
    from isograph_benchmark.real_data import locus_event_audit as lea

    monkeypatch.setattr(lea, "signal_root",
                        lambda arm="switch", max_snps=12000: Path("/nowhere") / f"g{max_snps}")
    # /nowhere has no parquet; the path named in the refusal is the point.
    with pytest.raises(SystemExit, match=r"/nowhere/g12000/all_introns/"):
        lea.load_nominations("susie", sqtl_arm=CONTEXT_DISTINCT_ARM)


def test_the_abf_layer_refuses_the_all_introns_arm():
    """coloc.abf is scored on the representative intron, so `--sqtl-arm all` there would
    be a claim about phenotypes that layer never fitted. Refuse rather than mis-tier."""
    from isograph_benchmark.real_data import locus_event_audit as lea

    with pytest.raises(SystemExit, match="no meaning"):
        lea.load_nominations("abf", sqtl_arm=CONTEXT_DISTINCT_ARM)


def test_a_raised_snp_guard_is_audited_from_its_own_root_never_the_primary(monkeypatch):
    """The MAX_SNPS=30000 aging__ad recovery is a scoped sensitivity arm. Its tiers (PICALM
    becomes signal-level there, and nowhere in the primary grid) must be read from its own
    results and written to their own directory, or a sensitivity result silently becomes
    the primary audit."""
    from isograph_benchmark.real_data import locus_event_audit as lea

    monkeypatch.setattr(lea, "signal_root",
                        lambda arm="switch", max_snps=12000: Path("/nowhere") / f"g{max_snps}")
    with pytest.raises(SystemExit, match=r"/nowhere/g30000/all_introns/"):
        lea.load_nominations("susie", sqtl_arm=CONTEXT_DISTINCT_ARM, max_snps=30000)

    monkeypatch.setattr(lea, "stage_out", lambda *a: Path("/nowhere").joinpath(*a[1:]))
    monkeypatch.setattr(lea, "ensure_dir", lambda p: p)
    assert lea.out_dir("susie", CONTEXT_DISTINCT_ARM).name == "susie_all_introns"
    assert (lea.out_dir("susie", CONTEXT_DISTINCT_ARM, 30000).name
            == "susie_all_introns__max_snps_30000")


def test_descriptors_join_on_the_whole_cell_so_a_gene_in_two_loci_is_not_duplicated(
        monkeypatch, tmp_path):
    """ZNF232 sits in two AD loci. Joining the hierarchy on (gene, trait, tissue) would
    attach each locus's estimator to both cells, doubling the cell and letting one locus's
    signal-level result certify the other."""
    from isograph_benchmark.real_data import locus_event_audit as lea

    cell = dict(analysis="aging__ad", trait="ad", gene="ENSG1", tissue="Brain_Cortex")
    pd.DataFrame([{**cell, "LOCUS_ID": "locus88_chr17", "symbol": "G", "PP4_sQTL": 0.9,
                   "PP4_eQTL": 0.1}]).to_parquet(tmp_path / "genes.parquet")
    pd.DataFrame([{**cell, "LOCUS_ID": "locus88_chr17", "PP4_sQTL": 0.9, "PP4_eQTL": 0.1},
                  {**cell, "LOCUS_ID": "locus239_chr17", "PP4_sQTL": 0.4, "PP4_eQTL": 0.1}]
                 ).to_parquet(tmp_path / "cells.parquet")
    pd.DataFrame([
        {**cell, "LOCUS_ID": "locus88_chr17", "modality": "sQTL", "estimator": "susie",
         "phenotype_id": "chr17:1:2:clu_1_-:ENSG1", "p12_min_call": 5e-6,
         "PP4_at_p12_sweep_min": 0.6, "fallback_reason": None},
        {**cell, "LOCUS_ID": "locus239_chr17", "modality": "sQTL", "estimator": "abf",
         "phenotype_id": None, "p12_min_call": np.nan, "PP4_at_p12_sweep_min": 0.2,
         "fallback_reason": "gwas_no_credible_set"},
        {**cell, "LOCUS_ID": "locus88_chr17", "modality": "eQTL", "estimator": "abf",
         "phenotype_id": None, "p12_min_call": np.nan, "PP4_at_p12_sweep_min": 0.0,
         "fallback_reason": "no_qtl_credible_set"},
    ]).to_parquet(tmp_path / "cells_hierarchy.parquet")
    monkeypatch.setattr(lea, "signal_root", lambda arm="switch", max_snps=12000: tmp_path)

    _, cells = lea.load_nominations("susie")
    assert len(cells) == 2
    by = cells.set_index("LOCUS_ID")
    assert by.at["locus88_chr17", "estimator_sQTL"] == "susie"
    assert by.at["locus239_chr17", "fallback_reason"] == "gwas_no_credible_set"
    assert by.at["locus88_chr17", "p12_min_call"] == pytest.approx(5e-6)


def test_a_hierarchy_written_before_the_descriptors_still_loads_without_inventing_them(
        monkeypatch, tmp_path):
    from isograph_benchmark.real_data import locus_event_audit as lea

    cell = dict(analysis="aging__lbd", trait="lbd", LOCUS_ID="locus08_chr4", gene="ENSG1",
                tissue="Brain_Cortex")
    pd.DataFrame([{**cell, "symbol": "G", "PP4_sQTL": 0.9, "PP4_eQTL": 0.1}]
                 ).to_parquet(tmp_path / "genes.parquet")
    pd.DataFrame([{**cell, "PP4_sQTL": 0.9, "PP4_eQTL": 0.1}]).to_parquet(tmp_path / "cells.parquet")
    pd.DataFrame([{**cell, "modality": "sQTL", "estimator": "susie", "phenotype_id": "x"}]
                 ).to_parquet(tmp_path / "cells_hierarchy.parquet")
    monkeypatch.setattr(lea, "signal_root", lambda arm="switch", max_snps=12000: tmp_path)

    _, cells = lea.load_nominations("susie")
    assert cells["estimator_sQTL"].iloc[0] == "susie"
    assert np.isnan(cells["p12_min_call"].iloc[0])
    assert cells["fallback_reason"].iloc[0] is None


def test_descriptors_describe_the_headline_tissue_not_a_more_flattering_one():
    """The quoted PP4 comes from `max_tissue_sQTL`. A second tissue that happens to be more
    prior-robust must not lend the headline its robustness; it is counted separately."""
    from isograph_benchmark.real_data.locus_event_audit import nomination_descriptors

    k = dict(analysis="aging__als", trait="als", LOCUS_ID="locus01_chr1", gene="ENSG1")
    nom = pd.DataFrame([{**k, "max_tissue_sQTL": "Brain_Cerebellum", "PP4_sQTL": 0.96}])
    cells = pd.DataFrame([
        {**k, "tissue": "Brain_Cerebellum", "PP4_sQTL": 0.96, "estimator_sQTL": "susie",
         "fallback_reason": None, "p12_min_call": 5e-6, "PP4_at_p12_sweep_min": 0.71},
        {**k, "tissue": "Brain_Cerebellar_Hemisphere", "PP4_sQTL": 0.94,
         "estimator_sQTL": "susie", "fallback_reason": None, "p12_min_call": 1e-6,
         "PP4_at_p12_sweep_min": 0.85},
        {**k, "tissue": "Brain_Cortex", "PP4_sQTL": 0.20, "estimator_sQTL": "abf",
         "fallback_reason": "no_qtl_credible_set", "p12_min_call": np.nan,
         "PP4_at_p12_sweep_min": 0.05},
    ])
    d = nomination_descriptors(nom, cells).iloc[0]
    assert d["headline_tissue"] == "Brain_Cerebellum"
    assert d["headline_estimator"] == "susie"
    assert d["prior_robustness"] == "intermediate"
    assert d["headline_PP4_at_p12_sweep_min"] == pytest.approx(0.71)
    assert d["n_tissue_sQTL_robust"] == 1


def test_an_abf_headline_carries_the_reason_susie_did_not_score_it():
    from isograph_benchmark.real_data.locus_event_audit import nomination_descriptors

    k = dict(analysis="aging__als", trait="als", LOCUS_ID="locus26_chr11", gene="ENSG1")
    nom = pd.DataFrame([{**k, "max_tissue_sQTL": "Brain_Cerebellum", "PP4_sQTL": 0.996}])
    cells = pd.DataFrame([{**k, "tissue": "Brain_Cerebellum", "PP4_sQTL": 0.996,
                           "estimator_sQTL": "abf", "fallback_reason": "gwas_no_credible_set",
                           "p12_min_call": 1e-6, "PP4_at_p12_sweep_min": 0.97}])
    d = nomination_descriptors(nom, cells).iloc[0]
    assert d["headline_fallback_reason"] == "gwas_no_credible_set"
    assert d["prior_robustness"] == "robust"


def test_the_abf_layer_refuses_a_snp_guard_it_never_applied():
    from isograph_benchmark.real_data import locus_event_audit as lea

    with pytest.raises(SystemExit, match="SNP guard"):
        lea.load_nominations("abf", max_snps=30000)


def test_an_abf_fallback_cell_may_not_certify_event_specificity():
    """The regression this gate was written for. PICALM's AD locus was dropped by the
    MAX_SNPS guard, so it has no coloc.susie row at all; its PP4 comes from the coloc.abf
    fallback, which is scored on GTEx's representative intron even inside the all-introns
    arm. Every other condition holds -- the curated interval IS a GTEx phenotype, the arm
    IS all-introns, a different intron IS named -- and the locus would still be tiered as
    event-specific off a test that never ran on any intron but one."""
    v = assign_tier(_introns((85974812, 85981129), chrom="chr11"),
                    _curated(gene="PICALM", trait="ad",
                             evidence_class="twas_association", chrom="chr11",
                             start=86026368, end=86031468),
                    curated_tested=True, sqtl_arm=CONTEXT_DISTINCT_ARM,
                    signal_level=False)
    assert v["tier"] == "disease_locus_splice_linked"
    assert v["signal_level"] is False


def test_a_tested_curated_event_that_stays_silent_is_context_distinct_not_merely_linked():
    """UNC13A. The cryptic-exon intron IS a GTEx phenotype, the all-introns arm fitted it,
    and it produced no colocalization -- while a different intron 9 kb upstream came back
    at PP4 ~0.996. Because the famous event was measured alongside the observed one and
    came back empty, the observed signal cannot be dismissed as a proxy for it. That is
    event specificity, and it is a stronger statement than `disease_locus_splice_linked`."""
    v = assign_tier(_introns(UNC13A_COLOC), _curated(), curated_tested=True,
                    sqtl_arm=CONTEXT_DISTINCT_ARM, signal_level=True)
    assert v["tier"] == "context_distinct_splice_colocalization"
    assert v["curated_event_tested"] is True
    assert v["curated_interval_matched"] is False


def test_an_untested_curated_event_stays_merely_linked():
    """The whole tier rests on the curated event having been on trial. If GTEx never
    carried a phenotype there, the non-match says nothing -- the event might colocalize
    beautifully and we would not know -- so no specificity may be claimed."""
    v = assign_tier(_introns(UNC13A_COLOC), _curated(), curated_tested=False,
                    sqtl_arm=CONTEXT_DISTINCT_ARM, signal_level=True)
    assert v["tier"] == "disease_locus_splice_linked"


def test_the_representative_arm_can_never_assign_context_distinct():
    """In the representative arm GTEx's grouped-permutation winner is the only intron
    fitted, so a non-representative curated event is never tested. Assigning the tier
    there would claim a negative result from a test that did not run -- exactly the error
    the all-introns arm was built to remove."""
    v = assign_tier(_introns(UNC13A_COLOC), _curated(), curated_tested=True,
                    sqtl_arm="representative", signal_level=True)
    assert v["tier"] == "disease_locus_splice_linked"


def test_context_distinct_never_outranks_an_actual_match():
    """If the curated event itself colocalized, that is recovery of the known mechanism
    and the new tier must not intercept it."""
    v = assign_tier(_introns(UNC13A_GTEX_CE), _curated(), curated_tested=True,
                    sqtl_arm=CONTEXT_DISTINCT_ARM, signal_level=True)
    assert v["tier"] == "known_mechanism_recovered"


def test_a_twas_capped_coordinate_match_is_not_context_distinct():
    """PICALM's TWAS anchor caps promotion, but a capped match still means the curated
    event DID colocalize. Deciding the new tier on `interval_matched` rather than on
    `matched` is what keeps those apart -- reading the evidence cap as 'the curated event
    was silent' would invert the finding."""
    cur = _curated(gene="PICALM", trait="ad", evidence_class="twas_association",
                   chrom="chr11", start=86026368, end=86031468)
    v = assign_tier(_introns((86026367, 86031469), chrom="chr11"), cur,
                    curated_tested=True, sqtl_arm=CONTEXT_DISTINCT_ARM,
                    signal_level=True)
    assert v["curated_interval_matched"] is True
    assert v["curated_matched"] is False
    assert v["tier"] == "disease_locus_splice_linked"


def test_context_distinct_needs_some_intron_to_have_colocalized():
    """A locus where nothing colocalized is `not_resolved`; the curated event being
    tested does not upgrade an empty result."""
    v = assign_tier([], _curated(), curated_tested=True,
                    sqtl_arm=CONTEXT_DISTINCT_ARM, signal_level=True)
    assert v["tier"] == "not_resolved"


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
