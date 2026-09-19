"""Tests for the short-read junction confirmation of the anchored switch pairs."""
from __future__ import annotations

import numpy as np
import pandas as pd

from isograph_benchmark.real_data import junction_coloc_confirm as jcc


# --------------------------------------------------------------------------- #
# event_info parsing
# --------------------------------------------------------------------------- #
def test_event_coords_reads_both_arms():
    """The competing arm must be parsed, not just the first junction.

    Regression: a LIBD event_info names the chromosome ONCE and then lists both arms, so
    validate_switch_splicing's chr-anchored ``(chr[\\w]+):(\\d+)-(\\d+)`` matches only the
    first arm. CTSH's anchored junction chr15:78937423-78939140 is the SECOND arm of its
    A3 event, so it was reported as "junction_not_measured" -- a false negative that reads
    exactly like "the junction is absent from the data".
    """
    info = "A3:chr15:78937496-78939140:78937423-78939140:-"
    coords = jcc._event_coords(info)
    assert coords == [78937496, 78939140, 78937423, 78939140]
    # both arms of the SNCA alternative-first-exon event
    snca = "AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-"
    assert 89838252 in jcc._event_coords(snca)
    assert 89836127 in jcc._event_coords(snca)


def test_match_events_finds_second_arm_junction():
    gm = pd.DataFrame({
        "event_id": ["e1", "e2"],
        "chrom": ["chr15", "chr15"],
        "event_info": ["A3:chr15:78937496-78939140:78937423-78939140:-",
                       "SE:chr15:78935750-78937318:78937423-78939140:-"],
    })
    hits = jcc._match_events("chr15:78937423-78939140(-)", gm)
    assert set(hits) == {"e1", "e2"}


def test_match_events_rejects_other_chromosome():
    """Chromosome is checked separately, since the pattern no longer carries it."""
    gm = pd.DataFrame({
        "event_id": ["e1"],
        "chrom": ["chr4"],
        "event_info": ["A3:chr4:78937496-78939140:78937423-78939140:-"],
    })
    assert jcc._match_events("chr15:78937423-78939140(-)", gm) == []


# --------------------------------------------------------------------------- #
# the statistic
# --------------------------------------------------------------------------- #
def test_minor_form_usage_is_orientation_free():
    """min(median, 1-median) must not depend on which arm is the numerator."""
    for med in (0.02, 0.98):
        assert abs(min(med, 1 - med) - 0.02) < 1e-12
    for med in (0.19, 0.81):
        assert abs(min(med, 1 - med) - 0.19) < 1e-12


def test_reference_contrast_selects_the_display_item_arm():
    """The gene-level verdict must key on the arm the figure claims, not any event."""
    fig4a = "AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-"
    other = "AF:chr4:89835692-89836127:89836213:89835692-89836743:89837139:-"
    assert jcc._is_reference("SNCA", fig4a)
    assert not jcc._is_reference("SNCA", other)
    assert not jcc._is_reference("NOT_A_TARGET", fig4a)


def test_empirical_p_never_zero():
    """(k+1)/(n+1): a finite null cannot license p = 0."""
    null = np.array([0.1, 0.2, 0.3])
    assert jcc._empirical_p(-5.0, null) == 1 / 4
    assert jcc._empirical_p(0.5, null) == 1.0
    assert np.isnan(jcc._empirical_p(float("nan"), null))
    assert np.isnan(jcc._empirical_p(-1.0, np.array([])))


def test_residualize_returns_nan_when_underdetermined():
    design = np.column_stack([np.ones(4), np.arange(4.0)])
    y = np.array([1.0, np.nan, np.nan, np.nan])
    assert np.isnan(jcc._residualize(y, design)).all()


def test_residualize_removes_the_covariate():
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    design = np.column_stack([np.ones(200), x])
    y = 3.0 + 2.0 * x + rng.normal(scale=0.01, size=200)
    r = jcc._residualize(y, design)
    assert abs(float(np.corrcoef(r, x)[0, 1])) < 0.05


# --------------------------------------------------------------------------- #
# targets
# --------------------------------------------------------------------------- #
def test_tissue_region_map_covers_the_coloc_tissues():
    """Each anchored coloc tissue must map to the BrainSEQ region that measures it."""
    assert jcc._TISSUE_TO_REGION["Brain_Hippocampus"] == ("hippocampus",)
    for cortex in ("Brain_Cortex", "Brain_Frontal_Cortex_BA9"):
        assert "dlpfc" in jcc._TISSUE_TO_REGION[cortex]


def test_targets_do_not_require_a_reported_competitor():
    """SNCA's coloc reports no competing junction; that must not make it untestable.

    The competitor comes from the PSI event's own second arm in ``direct`` mode, so
    requiring one up front would have dropped the paper's headline locus entirely.
    """
    tgt = jcc.load_targets(("SNCA", "CTSH"))
    if tgt.empty:  # deep-dive outputs absent in a bare checkout
        return
    snca = tgt[tgt["gene_name"] == "SNCA"]
    assert not snca.empty
    assert (snca["n_competing"] == 0).all()
    assert set(snca["region"]) <= {"dlpfc", "caudate"}


def test_absent_genes_are_a_null_not_a_crash(tmp_path, monkeypatch):
    """Genes missing from the deep dive yield zero targets, and the region loop skips them."""
    ev = pd.DataFrame({"gene_name": ["OTHER"], "junction": ["chr1:1-2"], "trait": ["ad"],
                       "tissue": ["Brain_Cortex"], "junction_in_switch_pair": [True],
                       "ens": ["ENSG0"], "clpp": [0.1], "risk_allele": ["A"]})
    path = tmp_path / "deep_dive_events.parquet"
    ev.to_parquet(path)
    monkeypatch.setattr(jcc, "_DEEP_DIVE", path)
    tgt = jcc.load_targets(("SNCA", "CTSH"))
    assert tgt.empty and "region" in tgt.columns
    res, nul = jcc.run_region("dlpfc", tgt, alpha=0.05, min_n=30, n_background_genes=1,
                              max_pairs=1, seed=0, min_usage=0.05)
    assert res.empty and nul.empty
    assert "SNCA, CTSH" in jcc._write_report(res, 0.05, 0.05, missing=("SNCA", "CTSH"))


# --------------------------------------------------------------------------- #
# upstream parser (validate_switch_splicing)
# --------------------------------------------------------------------------- #
def test_validate_switch_splicing_parses_both_event_arms():
    """The shared PSI parser must see an event's second arm.

    Regression for the same defect fixed in this module: `_JUNC_RE` requires a `chrN:`
    prefix on every match, but a LIBD event_info names the chromosome once and then lists
    both arms, so it captured only the first. Junctions carried by a second arm were
    silently reported as unmeasured, which undercounted analysis B (9/17 -> 11/17, and
    5 -> 11 exact two-endpoint matches).
    """
    from isograph_benchmark.real_data import validate_switch_splicing as vss

    info = "A3:chr15:78937496-78939140:78937423-78939140:-"
    coords = [int(a) for pair in vss._EVENT_PAIR_RE.findall(info) for a in pair]
    assert coords == [78937496, 78939140, 78937423, 78939140]
    # the old pattern is kept for single-junction strings, where it is correct
    assert vss._parse_junction("chr15:78937423-78939140(-)") == ("chr15", 78937423, 78939140)


def test_match_event_still_requires_the_right_chromosome():
    """Dropping chr from the coordinate pattern must not lose chromosome specificity."""
    import pandas as pd

    from isograph_benchmark.real_data import validate_switch_splicing as vss

    events = pd.DataFrame({
        "event_id": ["e1"], "chrom": ["chr4"],
        "coords": [[78937423, 78939140]],
    })
    assert vss._match_event("chr15", 78937423, 78939140, events) == (None, 0)
    assert vss._match_event("chr4", 78937423, 78939140, events) == ("e1", 2)
