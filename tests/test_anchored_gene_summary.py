"""The anchored-gene table leads with coloc, and never grades a switch on how rare it is."""
from __future__ import annotations

import numpy as np
import pandas as pd

from isograph_benchmark.real_data import anchored_gene_summary as ags


def _row(**kw):
    base = dict(
        gene="X", trait="ad", pp4_sqtl=0.95, pp4_eqtl=0.1,
        n_tissue_sqtl_coloc=5, n_tissue_eqtl_coloc=0,
        smr_sqtl_supported=4, smr_sqtl_tested=4, smr_sqtl_heidi_rejected=0,
        smr_eqtl_supported=0, smr_eqtl_tested=0, smr_eqtl_heidi_rejected=0,
        clpp=0.01,
        lr_confirmed=True, lr_expressed=True, lr_pairs_switchlike=2, lr_pairs_detected=4,
        ase_donors_both_forms=400, ase_usage_donors=498,
        ase_usage_min=0.045, ase_usage_median=0.2, ase_usage_max=0.39,
        ase_testable=3, ase_pairs=5, psi_verdict="junction not in catalogue",
    )
    base.update(kw)
    return pd.Series(base)


# --------------------------------------------------------------------------- #
# coloc is primary, CLPP is secondary
# --------------------------------------------------------------------------- #
def test_clpp_never_changes_the_genetics_call():
    """CLPP is noisy per locus and, for this set's two largest values, contradicted.

    TMEM175 carries CLPP 0.483 beside a coloc.abf PP4_sQTL of ~0; a table that ranked on
    CLPP would put it first. The call must be identical whatever CLPP says.
    """
    low = ags._genetics_call(_row(clpp=0.001))
    high = ags._genetics_call(_row(clpp=0.99))
    assert low == high == "splicing-specific"

    weak = _row(pp4_sqtl=0.0001, n_tissue_sqtl_coloc=0, smr_sqtl_supported=0,
                smr_sqtl_tested=0)
    assert ags._genetics_call(weak.copy()) == "weak coloc"
    weak_high_clpp = weak.copy()
    weak_high_clpp["clpp"] = 0.483            # the TMEM175 shape
    assert ags._genetics_call(weak_high_clpp) == "weak coloc"


def test_splicing_specific_requires_no_eqtl_instrument():
    """With an eQTL instrument also supported the gene is led, not specific."""
    assert ags._genetics_call(_row()) == "splicing-specific"
    shared = _row(smr_eqtl_tested=2, smr_eqtl_supported=2)
    assert ags._genetics_call(shared) == "splicing-led"


def test_tissue_counts_gate_the_splicing_contrast():
    """GTEx eQTL power exceeds sQTL power everywhere, so a lone PP4 is not a contrast."""
    r = _row(n_tissue_sqtl_coloc=1, n_tissue_eqtl_coloc=5)
    assert ags._genetics_call(r) == "colocalizes, not splicing-specific"


# --------------------------------------------------------------------------- #
# rarity is biology, not a defect
# --------------------------------------------------------------------------- #
def test_a_rare_but_well_detected_form_still_confirms():
    """BrainSEQ and GTEx are neurotypical tissue.

    A form that disease raises is expected to be a minority in controls, so gating on
    minor-form usage would discard exactly the disease-relevant switches. The call must
    key on whether both forms are observed, and in how many donors.
    """
    rare = _row(ase_usage_min=0.002, ase_usage_max=0.002, ase_donors_both_forms=495,
                ase_usage_donors=498, lr_confirmed=False, lr_pairs_switchlike=0)
    assert ags._orthogonal_call(rare).startswith("confirmed")

    # ...while a form seen in a handful of donors is not evidence of a switch
    thin = _row(ase_usage_min=0.48, ase_usage_max=0.48, ase_donors_both_forms=3,
                ase_usage_donors=4, lr_confirmed=False, lr_pairs_switchlike=0)
    assert ags._orthogonal_call(thin) == "tested, not confirmed"


def test_usage_threshold_is_reference_only():
    """The 0.05 figure exists for continuity with the PSI arm; nothing here applies it."""
    src = (ags.__file__.replace(".pyc", ".py"))
    with open(src) as fh:
        body = fh.read()
    assert "MIN_USAGE_REFERENCE" in body
    # it may be printed in the report, but must not appear in either call
    call_src = body[body.index("def _genetics_call"):body.index("# ---", body.index(
        "def _orthogonal_call"))]
    assert "MIN_USAGE_REFERENCE" not in call_src


def test_assays_are_counted_independently():
    none = _row(lr_confirmed=False, lr_pairs_switchlike=0, lr_expressed=False,
                ase_donors_both_forms=0, ase_testable=0, ase_pairs=0,
                psi_verdict="no BrainSEQ region")
    assert ags._orthogonal_call(none) == "not measurable"
    two = _row(psi_verdict="validated")          # long-read + ase + psi
    assert ags._orthogonal_call(two) == "confirmed (3 assays)"


# --------------------------------------------------------------------------- #
# reproducibility
# --------------------------------------------------------------------------- #
def test_provenance_names_every_input_and_the_commit():
    p = ags.provenance()
    assert p["git_commit"] and "generated_utc" in p and "git_dirty" in p
    assert p["generator"].endswith("anchored_gene_summary.py")
    paths = [f["path"] for f in p["inputs"]]
    for frag in ("coloc_isoform_events_combined", "smr_results", "gene_confirmation",
                 "junction_confirm", "pair_feasibility", "junction_allelic_counts"):
        assert any(frag in s for s in paths), frag
    for f in p["inputs"]:
        assert f.get("missing") or f["digest"].startswith("sha256")
    # thresholds travel with the numbers, so a reader can see what produced a call
    assert set(p["thresholds"]) >= {"PP4_STRONG", "MIN_DONORS_BOTH"}


def test_report_flags_a_dirty_tree():
    t = pd.DataFrame([_row()])
    t["tissue"] = "Cortex"
    t["go_invisible"] = True
    t["brainseq_switch_replicates"] = False
    t["genetics_call"] = "splicing-specific"
    t["orthogonal_call"] = "confirmed (2 assays)"
    t["ase_frac_donors_both"] = 0.8
    t["ase_usage_pairs"] = 3
    prov = dict(generated_utc="2026-09-20T00:00:00+00:00", git_commit="abcdef1234",
                git_dirty=True, pandas=pd.__version__, numpy=np.__version__,
                generator="x", thresholds={"PP4_STRONG": 0.8, "MIN_DONORS_BOTH": 30},
                inputs=[{"path": "a", "bytes": 1, "mtime": "t", "digest": "sha256:x"}])
    md = ags.report(t, prov)
    assert "not reproducible from that commit alone" in md
    assert "do not edit any number here by hand" in md.lower()
    assert "abcdef12" in md


# --------------------------------------------------------------------------- #
# the junction-recount columns must be about the anchored junction
# --------------------------------------------------------------------------- #
def test_pairs_carrying_matches_on_coordinates_with_convention_slack():
    """LeafCutter names the flanking exon bases; the recount stores the intron itself.

    PRDM2's colocalizing junction is chr1:13816570-13823159 in GTEx and appears as
    13816571-13823158 in the recount -- one base in at each end. Matching without slack
    silently finds nothing, which would read as "the junction is not measured".
    """
    jn = pd.DataFrame({
        "pair_id": ["A", "B", "C"],
        "chrom": ["chr1", "chr1", "chr2"],
        "start": [13816571, 13800000, 13816571],
        "end": [13823158, 13810000, 13823158],
    })
    want = {("chr1", 13816570, 13823159)}
    assert ags._pairs_carrying(jn, want) == {"A"}          # slack absorbs the offset
    assert ags._pairs_carrying(jn, {("chr1", 1, 2)}) == set()
    # a matching interval on another chromosome is not the same junction
    assert "C" not in ags._pairs_carrying(jn, want)


def test_anchored_junctions_parses_the_stranded_form():
    c = pd.DataFrame({
        "ens": ["E1", "E1", "E2"],
        "junction": ["chr1:13816570-13823159(+)", "chr1:1-2(+)", None],
    })
    got = ags._anchored_junctions(c)
    assert got["E1"] == {("chr1", 13816570, 13823159), ("chr1", 1, 2)}
    assert "E2" not in got       # a gene with no parseable junction is absent, not empty


def test_usage_is_a_range_because_the_partner_decides_the_ratio():
    """PRDM2's junction sits in 12 switch pairs spanning 0.045 to 0.39 minor-form usage.

    Electing one of them to stand for "the" usage hid the panel's own contrast (0.045)
    behind an unrelated pair (0.39), so the table reports the spread instead.
    """
    t = pd.DataFrame([_row()])
    for c in ("tissue", "go_invisible", "brainseq_switch_replicates", "genetics_call",
              "orthogonal_call"):
        t[c] = "x" if c in ("tissue", "genetics_call", "orthogonal_call") else False
    t["ase_frac_donors_both"] = 0.8
    t["ase_usage_pairs"] = 12
    md = ags.report(t)
    assert "0.045-0.390" in md
    assert "range" in md.lower()
