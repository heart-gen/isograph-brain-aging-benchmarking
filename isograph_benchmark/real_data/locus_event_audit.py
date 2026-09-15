"""Event-level audit: does a colocalizing locus name the splice event it claims to?

WHY THIS EXISTS
---------------
A colocalization posterior says a disease signal and a splicing signal share a causal
variant at a gene. It does NOT say which splice event, and the difference decides how the
result may be described. "UNC13A colocalizes with a brain sQTL" and "we recovered the
TDP-43-dependent cryptic-exon mechanism" look identical in a results table and are not the
same claim: the second is only true if the colocalizing intron IS the cryptic-exon intron.

The audit walks the chain the claim actually needs, per locus:

    GWAS credible signal
      -> GTEx sQTL credible signal
        -> the GTEx intron/junction phenotype that signal acts on
          -> the IsoGraph driver transcript(s) the intron belongs to

and then asks whether that intron matches a curated literature event
(`configs/known_splice_events.yaml`).

THE THING THIS WAS BUILT TO CATCH
---------------------------------
GTEx's sQTL phenotype for a gene is ONE grouped-permutation representative intron, chosen
from splicing data alone. That choice is not circular, but it is also not aimed at the
event anyone cares about: for UNC13A the representative intron in the tissues carrying the
colocalization is not the cryptic-exon intron at all. A gene-level posterior therefore
cannot, by construction, tell you whether the disease-relevant event was recovered. Only
running every intron of the gene and reporting which ones colocalize can.

THE TIER RULE, AND WHY `status` GATES IT
----------------------------------------
    known_mechanism_recovered   all four links converge on a REVIEWED curated event
    context_distinct_splice_colocalization
                                coloc holds at an event that is NOT the curated one, and
                                the curated event was itself TESTED and did not colocalize
    disease_locus_splice_linked coloc holds and an event is named, but not that event,
                                and the curated event was not shown to have been tested
    novel_splice_led_candidate  coloc holds, event named, no curated event at the locus
    not_resolved                coloc does not hold, or no intron could be named

`context_distinct_splice_colocalization` is the tier that needs justifying, because it is
a stronger claim than `disease_locus_splice_linked` and was added after the all-introns
arm ran (see AMENDMENT below). Three conditions guard it, and each closes a different way
of claiming a test that did not happen: the all-introns arm must be the one that produced
the numbers; GTEx must carry an sQTL phenotype at the curated intron; and the
colocalization must come from `coloc.susie` rather than the `coloc.abf` fallback, which is
scored on the representative intron whichever arm it sits in. The core of it is one fact:
whether the curated event was on trial. If GTEx never carried an sQTL phenotype at the curated intron, a non-match is
uninformative -- the event might colocalize beautifully and we would not know, so
`disease_locus_splice_linked` is all that is licensed. If GTEx DID carry that phenotype,
the arm DID fit it, and it still did not colocalize while a different intron of the same
gene did, then the negative is informative and event specificity is established: the
observed signal cannot be dismissed as a proxy for the famous event, because the famous
event was measured alongside it and came back empty. That is a different scientific
statement and it gets a different name.

The tier is reachable ONLY from the all-introns arm, and within it only from
signal-level cells. In the representative arm GTEx's grouped-permutation winner is the
single intron tested, so a curated non-representative event is never on trial and the
distinguishing fact is unavailable by construction. The `coloc.abf` fallback has the same
problem for the same reason -- it scores the representative intron alone -- and it fills
in wherever SuSiE found nothing, so it is the quiet path by which a locus that was never
fitted at all can arrive carrying a posterior. PICALM is exactly that case: its AD locus
exceeded the MAX_SNPS guard in 03c.coloc_gwas_susie.R and was dropped before any GWAS
SuSiE fit, so its PP4 is an abf number and no signal-level test has ever run there.

AMENDMENT, 2026-09-10
---------------------
The three-tier scheme above was pre-registered before the audit ran, and this tier was
added afterwards. Recorded plainly rather than folded in silently. What it is NOT is a
threshold moved to rescue a result: PP4_CALL, the interval tolerance, the `status` gate
and the `evidence_class` gate are all unchanged, and no locus that failed a gate now
passes one. What it IS: the original taxonomy could not express a distinction its own
motivating example turns out to need. `disease_locus_splice_linked` conflates "we found a
different event" with "we found a different event AND the curated one was tested and is
absent", and at UNC13A those are different findings. The tier splits that cell; it moves
no boundary. `known_mechanism_recovered` is untouched, so nothing can reach the top tier
that could not before.

Only a `reviewed` curated event -- coordinates independently verified in this repo against
GENCODE v47 -- can promote a locus to `known_mechanism_recovered`. A `proposed` record is
reported beside the verdict and can never change it. That is deliberate: the difference
between "the literature says this locus is spliced" and "we verified the coordinates of
the event and matched them" is exactly the difference between a positive control and a
story, and it should not rest on how confidently a note was written.

STAGES
------
  --stage audit   Walk the chain for the nominated loci and assign tiers.

The nomination source is switchable so the SAME audit runs unchanged against the current
`coloc.abf` grid and, later, against the signal-level `coloc.susie` nominations:
  --nominations abf    05_genetic_anchoring/_m/coloc_modality_contrast/{genes,cells}.parquet
  --nominations susie  05_genetic_anchoring/_m/coloc_signal_susie/{genes,cells}.parquet

Outputs under 05_genetic_anchoring/_m/locus_event_audit/<nominations>/.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from isograph_benchmark.paths import ensure_dir, rel, stage_out
from isograph_benchmark.real_data.coloc_modality_contrast import GTEX_BRAIN
from isograph_benchmark.real_data.coloc_signal_susie import (
    MAX_SNPS_PRIMARY,
    PRIOR_ROBUSTNESS,
    prior_robustness,
    results_dir,
    signal_root,
)

# Per-cell evidence-strength descriptors written by `coloc_signal_susie --stage meta`.
# Reported beside every nomination; none of them enters `assign_tier`.
_DESCRIPTORS = ("p12_min_call", "PP4_at_p12_sweep_min", "fallback_reason")

_EVENTS_YAML = rel("configs", "known_splice_events.yaml")
_GTF_CACHE = stage_out("tmp", "gencode.v47.primary_assembly.annotation.gtf_cache.parquet")

# Pre-specified before the audit was run.
PP4_CALL = 0.80
TIERS = ("known_mechanism_recovered", "context_distinct_splice_colocalization",
         "disease_locus_splice_linked", "novel_splice_led_candidate", "not_resolved")
# Two INDEPENDENT gates on the top tier, because they answer different questions.
# `status` asks whether the coordinates were verified here; `evidence_class` asks how
# strong the underlying literature claim is. A locus must clear both.
PROMOTING_STATUS = ("reviewed",)
# TWAS is deliberately absent. A TWAS tests association between genetically predicted
# splicing and the trait; it does NOT establish a shared causal variant, and LD-driven
# co-regulation at a locus produces the same signal. Colocalization is the more
# conservative test and is the bar this project holds itself to elsewhere, so a TWAS-only
# anchor caps a locus at `disease_locus_splice_linked` however well its coordinates match.
PROMOTING_EVIDENCE = ("functional_validation", "coloc_association")

# Only the all-introns arm puts a non-representative curated event on trial, so only it
# can distinguish "the curated event did not colocalize" from "the curated event was
# never tested". `context_distinct_splice_colocalization` rests entirely on that
# distinction and is therefore unreachable from any other arm.
CONTEXT_DISTINCT_ARM = "all"

GTEX_V11 = Path("/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11")
_SQTL_ALLPAIRS = "GTEx_Analysis_v11_sQTL_all_associations"


def out_dir(nominations: str, sqtl_arm: str = "representative",
            max_snps: int = MAX_SNPS_PRIMARY) -> Path:
    """One directory per (nomination source, sQTL arm, GWAS SNP guard). The arm changes both
    the numbers and the tiers, so writing both into `susie/` would leave whichever ran last;
    a non-default guard is a sensitivity arm and must not replace the primary audit."""
    leaf = nominations if sqtl_arm != CONTEXT_DISTINCT_ARM else f"{nominations}_all_introns"
    if int(max_snps) != MAX_SNPS_PRIMARY:
        leaf = f"{leaf}__max_snps_{int(max_snps)}"
    return ensure_dir(stage_out("anchoring", "locus_event_audit", leaf))


# --------------------------------------------------------------------------- #
# Curated events
# --------------------------------------------------------------------------- #
def load_curated_events() -> pd.DataFrame:
    d = yaml.safe_load(_EVENTS_YAML.read_text())
    ev = pd.DataFrame(d["events"])
    bad = set(ev["status"]) - {"reviewed", "proposed"}
    if bad:
        raise SystemExit(f"unknown status values in {_EVENTS_YAML}: {sorted(bad)}")
    if "evidence_class" not in ev.columns:
        raise SystemExit(f"every event needs an `evidence_class` in {_EVENTS_YAML}")
    known_ev = {"functional_validation", "coloc_association", "twas_association"}
    bad_ev = set(ev["evidence_class"].dropna()) - known_ev
    if bad_ev:
        raise SystemExit(f"unknown evidence_class values: {sorted(bad_ev)}")
    # A reviewed record without coordinates cannot gate anything, and silently letting it
    # through would make the tier depend on a record that names no interval.
    rev = ev[ev["status"] == "reviewed"]
    if rev[["chrom", "start", "end"]].isna().any().any():
        raise SystemExit(
            "a `reviewed` curated event is missing coordinates; either supply them or "
            "downgrade it to `proposed`")
    return ev


# --------------------------------------------------------------------------- #
# Link 3->4: intron -> transcripts
# --------------------------------------------------------------------------- #
def _parse_gtex_intron(phenotype_id: str) -> tuple[str, int, int] | None:
    """`chr19:17630750:17632782:clu_28403_-:ENSG...` -> (chrom, start, end)."""
    if not isinstance(phenotype_id, str):
        return None
    f = phenotype_id.split(":")
    if len(f) < 3:
        return None
    try:
        return f[0], int(f[1]), int(f[2])
    except ValueError:
        return None


def transcripts_spliced_by(chrom: str, start: int, end: int, gene_bare: str,
                           gtf: pd.DataFrame, tol: int = 2) -> list[str]:
    """Transcripts of a gene whose spliced intron matches (chrom, start, end).

    LeafCutter intron bounds and GENCODE exon boundaries can differ by a base, so the
    match allows a small tolerance -- the same convention `coloc_isoform_events` uses.
    """
    g = gtf[(gtf["gene_bare"] == gene_bare) & (gtf["feature"] == "exon")
            & (gtf["chrom"] == chrom)]
    hits = []
    for tx, ex in g.groupby("transcript_id", sort=True):
        e = ex.sort_values("start")
        istart = e["end"].to_numpy()[:-1] + 1
        iend = e["start"].to_numpy()[1:] - 1
        if np.any((np.abs(istart - start) <= tol) & (np.abs(iend - end) <= tol)):
            hits.append(str(tx))
    return hits


def _load_gtf() -> pd.DataFrame:
    if not _GTF_CACHE.exists():
        raise SystemExit(f"missing GENCODE v47 cache: {_GTF_CACHE}")
    g = pd.read_parquet(_GTF_CACHE, columns=["transcript_id", "gene_id", "chrom",
                                             "strand", "feature", "start", "end"])
    g = g[g["feature"] == "exon"].copy()
    g["gene_bare"] = g["gene_id"].astype(str).str.split(".").str[0]
    return g


# --------------------------------------------------------------------------- #
# Was the curated event actually on trial?
# --------------------------------------------------------------------------- #
def curated_event_testability(curated: pd.DataFrame, tissues=None,
                              cache: Path | None = None) -> pd.DataFrame:
    """Per curated event: was a GTEx sQTL phenotype carried at that exact intron, and did
    GTEx fine-map a credible set there?

    These are two different facts and the tier depends on the first, not the second.

    `tested` -- a phenotype with those bounds exists in the GTEx sQTL all-pairs for at
    least one brain tissue. This is what licenses reading a non-match as a result: the
    intron was measured and entered the fit. Read from the all-pairs release rather than
    from our own output, because our output only contains introns that produced a signal
    pair -- inferring "tested" from it would be circular, and would silently call every
    event that failed to colocalize "untested".

    `gtex_credible_set` -- GTEx's own in-sample-LD SuSiE found a credible set for that
    phenotype. Reported, never gating. Its absence says there is no fine-mappable sQTL at
    the event in healthy brain, which for a pathology-dependent cryptic exon is the
    expected biology rather than a failure; but a phenotype can be tested and colocalize
    on our side without GTEx having called a credible set, so gating on it would discard
    real results.
    """
    tissues = tuple(tissues) if tissues is not None else GTEX_BRAIN
    if cache is not None and cache.exists():
        return pd.read_parquet(cache)

    import pyarrow.compute as pc
    import pyarrow.dataset as pds

    gtex_cs_f = stage_out("anchoring.coloc_signal") / "gtex_credible_sets.parquet"
    gcs = (pd.read_parquet(gtex_cs_f, columns=["phenotype_id", "tissue", "kind"])
           if gtex_cs_f.exists() else
           pd.DataFrame(columns=["phenotype_id", "tissue", "kind"]))
    gcs = gcs[gcs["kind"] == "sQTL"] if len(gcs) else gcs

    ev = curated[curated[["chrom", "start", "end"]].notna().all(axis=1)].copy()
    if ev.empty:
        return pd.DataFrame()
    ev["gene_bare"] = ev["gene_id"].astype(str).str.split(".").str[0]

    # Each all-pairs file is scanned ONCE for every curated event on its chromosome. The
    # filter is a substring match on phenotype_id, which arrow cannot satisfy from
    # statistics, so the phenotype column of a multi-GB file is read end to end -- doing
    # that per event rather than per file would re-read the same column N times.
    found: dict[tuple, set] = {}
    for chrom, grp in ev.groupby("chrom"):
        chrnum = str(chrom).removeprefix("chr")
        bares = sorted(set(grp["gene_bare"]))
        expr = pc.match_substring(pds.field("phenotype_id"), bares[0])
        for b in bares[1:]:
            expr = expr | pc.match_substring(pds.field("phenotype_id"), b)
        for tis in tissues:
            f = GTEX_V11 / _SQTL_ALLPAIRS / f"{tis}.v11.cis_sqtl.allpairs.chr{chrnum}.parquet"
            if not f.exists():
                continue
            t = pds.dataset(str(f)).to_table(columns=["phenotype_id"], filter=expr)
            found[(str(chrom), tis)] = set(t.column("phenotype_id").to_pylist())

    rows = []
    for cr in ev.itertuples(index=False):
        chrom, cstart, cend = str(cr.chrom), int(cr.start), int(cr.end)
        _tol = getattr(cr, "match_tolerance", 2)
        tol = 2 if (_tol is None or pd.isna(_tol)) else int(_tol)
        tested_in, cs_in, pids = [], [], set()
        for tis in tissues:
            hit = False
            for pid in found.get((chrom, tis), ()):
                parsed = _parse_gtex_intron(pid)
                if parsed is None or parsed[0] != chrom:
                    continue
                # LeafCutter bounds sit 1 bp outside the GENCODE intron at each end, the
                # same convention `interval_matches` is built around; reuse it so the
                # testability probe and the tier cannot disagree about what "this intron"
                # means.
                if interval_matches(parsed[1], parsed[2], cstart, cend,
                                    mode="exact", tol=tol):
                    hit = True
                    pids.add(pid)
            if hit:
                tested_in.append(tis)
                if len(gcs) and ((gcs["tissue"] == tis)
                                 & gcs["phenotype_id"].isin(pids)).any():
                    cs_in.append(tis)
        rows.append({
            "gene": str(cr.gene), "trait": str(cr.trait),
            "curated_interval": f"{chrom}:{cstart}-{cend}",
            "tested": bool(tested_in),
            "n_tissues_tested": len(tested_in),
            "tissues_tested": ";".join(tested_in),
            "gtex_credible_set": bool(cs_in),
            "n_tissues_gtex_cs": len(cs_in),
            "tissues_gtex_cs": ";".join(cs_in),
            "gtex_phenotype_ids": ";".join(sorted(pids)),
        })
    out = pd.DataFrame(rows)
    if cache is not None and len(out):
        out.to_parquet(cache, index=False)
    return out


# --------------------------------------------------------------------------- #
# Nominations
# --------------------------------------------------------------------------- #
def load_nominations(source: str, call: float = PP4_CALL,
                     sqtl_arm: str = "representative",
                     max_snps: int = MAX_SNPS_PRIMARY) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(gene-level nominations, per-tissue cells) for a nomination source.

    The two sQTL arms write the same file names, so the all-introns results live in their
    own subdirectory rather than on top of the representative ones. Resolving that here
    is what stops `--sqtl-arm all` from tiering the representative nominations: the arm
    would then govern the tier while the numbers came from the other arm, which is the
    one way this audit could silently certify event specificity it never tested.

    `max_snps` selects the GWAS SNP-guard root the same way: a raised guard is a scoped
    sensitivity arm whose results live beside, never on top of, the primary grid.
    """
    if source == "abf":
        base = stage_out("anchoring.coloc_modality")
        if sqtl_arm == CONTEXT_DISTINCT_ARM:
            raise SystemExit(
                "--sqtl-arm all has no meaning for --nominations abf: the coloc.abf layer "
                "is scored on the representative intron only.")
        if int(max_snps) != MAX_SNPS_PRIMARY:
            raise SystemExit(
                "--max-snps has no meaning for --nominations abf: the coloc.abf layer "
                "applies no GWAS SuSiE SNP guard.")
    elif source == "susie":
        base = results_dir(signal_root("switch", max_snps),
                           "all" if sqtl_arm == CONTEXT_DISTINCT_ARM else "representative")
    else:
        raise SystemExit(f"unknown nomination source {source!r}")
    gf, cf = base / "genes.parquet", base / "cells.parquet"
    for f in (gf, cf):
        if not f.exists():
            raise SystemExit(f"missing {f}; run that layer's --stage meta first")
    genes = pd.read_parquet(gf)
    cells = pd.read_parquet(cf)

    # Which estimator actually produced each cell's sQTL posterior. `cells.parquet` keeps
    # only the posteriors; the provenance lives beside it in `cells_hierarchy.parquet`.
    # Carrying it is not bookkeeping: a `coloc.abf` fallback row is scored on GTEx's
    # representative intron ALONE even inside the all-introns arm, so it never put a
    # non-representative curated event on trial. Tiering such a row as event-specific
    # would credit it with a test that did not run -- and because the fallback silently
    # fills the cells where SuSiE found nothing, that error is invisible downstream.
    #
    # The evidence-strength descriptors ride along from the same file: the prior range a
    # call survives, and why an abf cell was not scored by SuSiE. A table written before
    # those columns existed still loads; the descriptors are then absent, never invented.
    hf = base / "cells_hierarchy.parquet"
    if hf.exists():
        h = pd.read_parquet(hf)
        h = (h[h["modality"] == "sQTL"]
             .rename(columns={"estimator": "estimator_sQTL"}))
        # On the whole cell identity, not (gene, trait, tissue): a gene can sit in two loci
        # for one trait (ZNF232/AD, locus88 and locus239), and the narrower join duplicates
        # the cell and every per-tissue count built on it.
        on = [c for c in ("analysis", "trait", "LOCUS_ID", "gene", "tissue")
              if c in cells.columns and c in h.columns]
        keep = [*on, "estimator_sQTL", "phenotype_id",
                *(c for c in _DESCRIPTORS if c in h.columns)]
        cells = cells.merge(h[keep], on=on, how="left")
    else:
        cells["estimator_sQTL"] = None
    for c in _DESCRIPTORS:
        if c not in cells.columns:
            cells[c] = None if c == "fallback_reason" else np.nan

    nom = genes[genes["PP4_sQTL"] >= call].copy()
    return nom, cells


def _sqtl_phenotype_for(gene: str, tissue: str, srep: pd.DataFrame) -> str | None:
    m = srep[(srep["gene"] == gene) & (srep["tissue"] == tissue)]
    return None if m.empty else str(m["phenotype_id"].iloc[0])


# --------------------------------------------------------------------------- #
# The audit
# --------------------------------------------------------------------------- #
def interval_matches(istart: int, iend: int, cstart: int, cend: int,
                     mode: str = "exact", tol: int = 2) -> bool:
    """Does an observed intron correspond to a curated interval?

    `exact` requires BOTH endpoints within tolerance. The tolerance is not slack: it is
    the 1 bp-per-end convention difference between LeafCutter intron bounds and GENCODE
    exon boundaries. `overlap` is far looser and is only correct when the curated
    interval is a REGION rather than a specific intron -- introns within one LeafCutter
    cluster are alternative splice choices that routinely share an endpoint, so an
    overlap rule would let a different splice event of the same gene certify as the
    curated one. That is the error this audit exists to prevent.
    """
    if mode == "exact":
        return abs(istart - cstart) <= tol and abs(iend - cend) <= tol
    return istart <= cend + tol and iend >= cstart - tol


def assign_tier(introns: list, cur: pd.DataFrame, curated_tested: bool = False,
                sqtl_arm: str = "representative",
                signal_level: bool = False) -> dict:
    """Tier a locus from its colocalizing introns and the curated events at that locus.

    Two INDEPENDENT gates guard the top tier and a locus must clear both:
      * `status` must be `reviewed`  -- the coordinates were verified in this repo.
      * `evidence_class` must be functional or coloc-based -- a TWAS association does not
        establish a shared causal variant, so it may never certify a mechanism.
    `curated_interval_matched` is reported separately so that a coordinate match capped by
    the evidence gate stays visible instead of looking like no match at all.

    `curated_tested`, `sqtl_arm` and `signal_level` govern only the middle tier. A locus
    whose curated event was on trial and did not colocalize, while another intron of the
    same gene did, earns `context_distinct_splice_colocalization` -- the negative is
    informative, so the observed signal is demonstrably not a proxy for the curated event.
    All three conditions are required, and each rules out a different way of claiming a
    test that did not happen:
      * `sqtl_arm` must be the all-introns arm -- elsewhere only GTEx's permutation winner
        is fitted, so a non-representative curated event is never on trial.
      * `curated_tested` -- GTEx must carry an sQTL phenotype at the curated intron at all.
      * `signal_level` -- the colocalization must come from `coloc.susie`, not from the
        `coloc.abf` fallback. The fallback is scored on the representative intron even
        inside the all-introns arm, so it tests one intron no matter which arm it sits in.
    """
    matched = interval_matched = False
    matched_status = matched_evidence = None
    if len(cur) and introns:
        for cr in cur.itertuples(index=False):
            # A record with no interval (a `proposed` one that was never pinned) cannot
            # be matched against. pandas renders a missing YAML value as NaN, not None.
            if pd.isna(cr.chrom) or pd.isna(cr.start) or pd.isna(cr.end):
                continue
            _tol = getattr(cr, "match_tolerance", 2)
            tol = 2 if (_tol is None or pd.isna(_tol)) else int(_tol)
            _mode = getattr(cr, "match_mode", "exact")
            mode = "exact" if (_mode is None or pd.isna(_mode)) else str(_mode)
            for chrom, istart, iend, _pid in introns:
                if chrom != cr.chrom:
                    continue
                if not interval_matches(istart, iend, int(cr.start), int(cr.end),
                                        mode=mode, tol=tol):
                    continue
                interval_matched = True
                matched_status = str(cr.status)
                matched_evidence = str(getattr(cr, "evidence_class", ""))
                if (matched_status in PROMOTING_STATUS
                        and matched_evidence in PROMOTING_EVIDENCE):
                    matched = True

    # The curated event was on trial and came back empty, while some OTHER intron of the
    # gene colocalized. Note this is decided on `interval_matched`, not on `matched`: a
    # coordinate match capped by the evidence gate (PICALM's TWAS anchor) still means the
    # curated event DID colocalize, which is not this situation at all.
    context_distinct = bool(
        introns and len(cur) and not interval_matched and curated_tested
        and sqtl_arm == CONTEXT_DISTINCT_ARM and signal_level)

    if not introns:
        tier = "not_resolved"
    elif matched:
        tier = "known_mechanism_recovered"
    elif context_distinct:
        tier = "context_distinct_splice_colocalization"
    elif len(cur):
        tier = "disease_locus_splice_linked"
    else:
        tier = "novel_splice_led_candidate"
    return {"tier": tier, "curated_matched": matched,
            "curated_interval_matched": interval_matched,
            "curated_status": matched_status,
            "curated_evidence_class": matched_evidence,
            "curated_event_tested": bool(curated_tested),
            "signal_level": bool(signal_level)}


def nomination_descriptors(nom: pd.DataFrame, cells: pd.DataFrame,
                           call: float = PP4_CALL) -> pd.DataFrame:
    """Evidence strength of each nominated gene x trait, read at its HEADLINE tissue.

    The headline tissue is the one the gene-level PP4_sQTL was taken from
    (`max_tissue_sQTL`), so the descriptors describe the number that is actually quoted,
    not a more flattering tissue. `n_tissue_sQTL_robust` restates the tissue-consistency
    count at the most conservative prior: tissues whose call still holds at p12 = 1e-6.
    None of this enters `assign_tier`.
    """
    keys = [k for k in ("analysis", "trait", "LOCUS_ID", "gene")
            if k in nom.columns and k in cells.columns]
    cols = [*keys, "headline_tissue", "headline_estimator", "headline_fallback_reason",
            "headline_p12_min_call", "headline_PP4_at_p12_sweep_min", "prior_robustness",
            "n_tissue_sQTL_robust"]
    rows = []
    for r in nom.itertuples(index=False):
        m = np.ones(len(cells), dtype=bool)
        for k in keys:
            m &= (cells[k] == getattr(r, k)).to_numpy()
        c = cells[m]
        tis = getattr(r, "max_tissue_sQTL", None)
        head = c[c["tissue"] == tis]
        h = head.iloc[0] if len(head) else None
        p12c = h["p12_min_call"] if h is not None else np.nan
        p12c = float(p12c) if pd.notna(p12c) else np.nan
        called = c[c["PP4_sQTL"] >= call]
        rows.append({
            **{k: getattr(r, k) for k in keys},
            "headline_tissue": tis,
            "headline_estimator": None if h is None else h["estimator_sQTL"],
            "headline_fallback_reason": None if h is None else h["fallback_reason"],
            "headline_p12_min_call": p12c,
            "headline_PP4_at_p12_sweep_min": (np.nan if h is None
                                              else h["PP4_at_p12_sweep_min"]),
            "prior_robustness": prior_robustness(p12c),
            "n_tissue_sQTL_robust": int(sum(
                prior_robustness(float(x)) == "robust"
                for x in called["p12_min_call"] if pd.notna(x))),
        })
    return pd.DataFrame(rows, columns=cols)


def run_audit(nominations: str = "abf", call: float = PP4_CALL,
              sqtl_arm: str = "representative",
              max_snps: int = MAX_SNPS_PRIMARY) -> Path:
    dest = out_dir(nominations, sqtl_arm, max_snps)
    curated = load_curated_events()
    nom, cells = load_nominations(nominations, call=call, sqtl_arm=sqtl_arm,
                                  max_snps=max_snps)
    gtf = _load_gtf()

    # Which curated events were actually on trial. Cached per arm: the probe reads the
    # GTEx all-pairs, which does not change, but the tier it feeds does depend on the arm.
    test_f = dest / "curated_event_testability.parquet"
    testability = curated_event_testability(curated, cache=test_f)
    tested_keys = set()
    if len(testability):
        tested_keys = {(r.gene, r.trait) for r in
                       testability[testability["tested"]].itertuples(index=False)}

    srep_f = stage_out("anchoring.coloc_modality") / "sqtl_representative.parquet"
    srep = pd.read_parquet(srep_f) if srep_f.exists() else pd.DataFrame(
        columns=["gene", "tissue", "phenotype_id"])

    print(f"  nominations ({nominations}): {len(nom):,} gene x trait cells at "
          f"PP4_sQTL >= {call}")
    print(f"  curated events: "
          f"{(curated['status'] == 'reviewed').sum()} reviewed, "
          f"{(curated['status'] == 'proposed').sum()} proposed")

    rows = []
    for r in nom.itertuples(index=False):
        gene = str(r.gene)
        sym = str(getattr(r, "symbol", "") or gene)
        trait = str(r.trait)
        # Link 2->3: the tissues where this gene's sQTL actually colocalizes, and the
        # intron phenotype each of those used.
        m = (cells["gene"] == gene) & (cells["trait"] == trait) & (cells["PP4_sQTL"] >= call)
        # And the locus: a gene nominated at one locus must not borrow tissues from another.
        if "LOCUS_ID" in cells.columns and pd.notna(getattr(r, "LOCUS_ID", None)):
            m &= cells["LOCUS_ID"] == r.LOCUS_ID
        c = cells[m]
        introns, tissues, tx_all = [], [], set()
        n_signal_level = 0
        for cr in c.itertuples(index=False):
            tis = str(cr.tissue)
            pid = (str(cr.phenotype_id) if "phenotype_id" in c.columns
                   and pd.notna(getattr(cr, "phenotype_id", None))
                   else _sqtl_phenotype_for(gene, tis, srep))
            parsed = _parse_gtex_intron(pid) if pid else None
            if parsed is None:
                continue
            chrom, istart, iend = parsed
            # An abf-fallback cell still names an intron worth REPORTING -- it is GTEx's
            # representative one -- but it must not certify event specificity, so the
            # count is kept separate from the intron list rather than filtering it.
            if str(getattr(cr, "estimator_sQTL", "")) == "susie":
                n_signal_level += 1
            tissues.append(tis)
            introns.append((chrom, istart, iend, pid))
            tx_all.update(transcripts_spliced_by(chrom, istart, iend, gene, gtf))

        # Link 4 and the curated match. Extracted into `assign_tier` so the rule that
        # decides whether a locus may be called a recovered mechanism is testable on its
        # own, rather than only reachable by running the whole audit.
        cur = curated[(curated["gene"] == sym) & (curated["trait"] == trait)]
        verdict = assign_tier(introns, cur, curated_tested=(sym, trait) in tested_keys,
                              sqtl_arm=sqtl_arm, signal_level=n_signal_level > 0)
        matched = verdict["curated_matched"]
        interval_matched = verdict["curated_interval_matched"]
        matched_status = verdict["curated_status"]
        matched_evidence = verdict["curated_evidence_class"]
        tier = verdict["tier"]

        rows.append({
            "analysis": getattr(r, "analysis", None), "trait": trait,
            "gene": gene, "symbol": sym,
            "LOCUS_ID": getattr(r, "LOCUS_ID", None),
            "PP4_sQTL": float(r.PP4_sQTL), "PP4_eQTL": float(r.PP4_eQTL),
            "n_tissue_sQTL_coloc": int(getattr(r, "n_tissue_sQTL_coloc", len(c))),
            "n_tissue": int(getattr(r, "n_tissue", 0)),
            "n_introns_named": len(set((a, b, c_) for a, b, c_, _ in introns)),
            "colocalizing_introns": ";".join(sorted({f"{a}:{b}-{c_}"
                                                     for a, b, c_, _ in introns})),
            "tissues": ";".join(sorted(set(tissues))),
            "n_driver_transcripts": len(tx_all),
            "driver_transcripts": ";".join(sorted(tx_all)),
            "curated_event": bool(len(cur)),
            "curated_status": matched_status if matched_status else (
                str(cur["status"].iloc[0]) if len(cur) else None),
            "curated_evidence_class": (matched_evidence if matched_evidence else
                                       (str(cur["evidence_class"].iloc[0]) if len(cur) else None)),
            "curated_interval_matched": interval_matched,
            "curated_matched": matched,
            "curated_event_tested": verdict["curated_event_tested"],
            "n_cells_signal_level": n_signal_level,
            "estimator": "susie" if n_signal_level else "abf",
            "sqtl_arm": sqtl_arm,
            "tier": tier,
        })

    audit = pd.DataFrame(rows).sort_values(["tier", "PP4_sQTL"], ascending=[True, False])
    # Evidence strength beside every nomination: the prior range its headline call
    # survives, and why an abf headline was not scored by SuSiE. Descriptors, not gates.
    if len(audit):
        dk = [k for k in ("analysis", "trait", "LOCUS_ID", "gene")
              if k in nom.columns and k in cells.columns]
        audit = audit.merge(nomination_descriptors(nom, cells, call=call), on=dk, how="left")
    audit["max_snps"] = int(max_snps)
    audit.to_parquet(dest / "locus_event_audit.parquet", index=False)
    _write_report(dest, audit, curated, nominations, call, testability, sqtl_arm,
                  max_snps)
    print(f"\n  tiers: {audit['tier'].value_counts().to_dict()}")
    print(f"  wrote {dest}")
    return dest


def _write_report(dest: Path, audit: pd.DataFrame, curated: pd.DataFrame,
                  nominations: str, call: float,
                  testability: pd.DataFrame | None = None,
                  sqtl_arm: str = "representative",
                  max_snps: int = MAX_SNPS_PRIMARY) -> None:
    L: list[str] = []
    A = L.append
    A("# Locus event audit")
    A("")
    if int(max_snps) != MAX_SNPS_PRIMARY:
        A(f"**SENSITIVITY ARM: GWAS SuSiE SNP guard {int(max_snps):,}** (primary "
          f"{MAX_SNPS_PRIMARY:,}). A tier assigned here to a locus the primary grid could "
          "not fit is a sensitivity result, and may not enter the biological narrative "
          "without a locus-specific LD audit (`08d.locus_ld_robustness.sh`).")
        A("")
    A(f"Nomination source: `{nominations}` at `PP4_sQTL >= {call}`.")
    A("")
    A("A colocalization posterior names a gene, not an event. This audit walks "
      "GWAS signal -> sQTL signal -> intron phenotype -> driver transcript, and asks "
      "whether the colocalizing intron matches a curated literature event.")
    A("")
    A("## Curated events")
    A("")
    A("| gene | trait | status | evidence class | interval | can promote |")
    A("|---|---|---|---|---|---|")
    for r in curated.itertuples(index=False):
        iv = (f"{r.chrom}:{int(r.start):,}-{int(r.end):,}"
              if pd.notna(r.chrom) and pd.notna(r.start) and pd.notna(r.end)
              else "not pinned")
        ok = (str(r.status) in PROMOTING_STATUS
              and str(getattr(r, "evidence_class", "")) in PROMOTING_EVIDENCE)
        A(f"| {r.gene} | {r.trait} | **{r.status}** | `{r.evidence_class}` | {iv} | "
          f"{'yes' if ok else '**no**'} |")
    A("")
    A("Promotion to `known_mechanism_recovered` needs BOTH gates. `status: reviewed` "
      "means the coordinates were independently verified here against GENCODE v47. "
      "`evidence_class` must be `functional_validation` or `coloc_association`: a TWAS "
      "association does **not** establish a shared causal variant -- LD-driven "
      "co-regulation yields the same signal -- so a TWAS-only anchor caps a locus at "
      "`disease_locus_splice_linked` however well its coordinates match. Where that cap "
      "binds, `curated_interval_matched` still records whether the interval lined up.")
    A("")
    A("## Was the curated event on trial?")
    A("")
    A(f"sQTL arm: `{sqtl_arm}`.")
    A("")
    if testability is None or testability.empty:
        A("_Testability was not probed._")
    else:
        A("| gene | trait | curated interval | tested in GTEx | GTEx credible set |")
        A("|---|---|---|---|---|")
        for r in testability.itertuples(index=False):
            A(f"| {r.gene} | {r.trait} | `{r.curated_interval}` | "
              f"{'**yes** (' + str(r.n_tissues_tested) + '/13 tissues)' if r.tested else 'no'} | "
              f"{'yes (' + str(r.n_tissues_gtex_cs) + '/13)' if r.gtex_credible_set else '**no**'} |")
        A("")
        A("`tested` is read from the GTEx sQTL all-pairs release, not from this "
          "pipeline's output: our output only contains introns that produced a signal "
          "pair, so inferring testability from it would call every event that failed to "
          "colocalize \"untested\". It is what licenses reading a non-match as a result "
          "rather than as a gap, and it is the fact that separates "
          "`context_distinct_splice_colocalization` from `disease_locus_splice_linked`.")
        A("")
        A("`GTEx credible set` is reported and never gates. Its absence means there is no "
          "fine-mappable sQTL at that event in healthy brain -- which for an event that "
          "only appears under pathology is the expected biology, not a failed test.")
    A("")
    A("## Tier counts")
    A("")
    for t in TIERS:
        A(f"- `{t}`: {int((audit['tier'] == t).sum())}")
    A("")
    if sqtl_arm != CONTEXT_DISTINCT_ARM:
        A(f"`context_distinct_splice_colocalization` is unreachable in the `{sqtl_arm}` "
          "arm and its count above is 0 by construction, not by result: GTEx's "
          "grouped-permutation winner is the only intron tested, so a curated "
          "non-representative event is never on trial. Only "
          "`--sqtl-arm all` can assign it.")
        A("")
    if len(audit) and "prior_robustness" in audit.columns:
        A("## Evidence strength per nomination")
        A("")
        A("Two descriptors travel with every nominated PP4, read at the headline tissue the "
          "gene-level number was taken from. Neither gates a tier.")
        A("")
        A("- **Prior robustness**: the smallest `p12` in the sweep at which PP4 still clears "
          f"{call}. `robust` holds at 1e-6, `intermediate` from 5e-6, `primary_prior` only "
          "from the pre-specified 1e-5 upward. PP4 is monotone in p12, so this is the range "
          "the call survives rather than one arbitrary prior.")
        A("- **Why coloc.abf**, for a headline `coloc.susie` did not score: the first place "
          "the cell left the signal-level pipeline. `gwas_locus_over_max_snps` was never "
          "tested at signal level. `gwas_no_credible_set` was tested and the GWAS did not "
          "fine-map, which weakens any colocalization claimed there. "
          "`qtl_cs_not_matching_gtex` is the reference-LD artefact the agreement filter "
          "exists to remove.")
        A("")
        xt = pd.crosstab(audit["headline_estimator"].fillna("none"), audit["prior_robustness"])
        labels = [lb for lb in PRIOR_ROBUSTNESS if lb in xt.columns]
        A("| headline estimator | " + " | ".join(f"`{lb}`" for lb in labels) + " | total |")
        A("|---|" + "---|" * (len(labels) + 1))
        for est, row in xt.iterrows():
            A(f"| {est} | " + " | ".join(str(int(row[lb])) for lb in labels)
              + f" | {int(row.sum())} |")
        A("")
        fb = audit.loc[audit["headline_estimator"] == "abf",
                       "headline_fallback_reason"].value_counts()
        if len(fb):
            A("| why coloc.abf | nominations |")
            A("|---|---|")
            for k, v in fb.items():
                A(f"| `{k}` | {v} |")
            A("")
        A("| gene | trait | PP4 sQTL | headline tissue | estimator | why abf | "
          "PP4 at p12=1e-6 | calls from p12 | prior | robust tissues |")
        A("|---|---|---|---|---|---|---|---|---|---|")
        for r in audit.sort_values(["headline_estimator", "PP4_sQTL"],
                                   ascending=[False, False]).itertuples(index=False):
            why = r.headline_fallback_reason
            lo, p = r.headline_PP4_at_p12_sweep_min, r.headline_p12_min_call
            A(f"| {r.symbol} | {r.trait} | {r.PP4_sQTL:.3f} | {r.headline_tissue} | "
              f"{r.headline_estimator} | {f'`{why}`' if isinstance(why, str) else '—'} | "
              f"{f'{lo:.3f}' if pd.notna(lo) else '—'} | {f'{p:g}' if pd.notna(p) else '—'} | "
              f"`{r.prior_robustness}` | {r.n_tissue_sQTL_robust}/{r.n_tissue_sQTL_coloc} |")
        A("")
    A("## Loci with a curated event")
    A("")
    sub = audit[audit["curated_event"]]
    if sub.empty:
        A("_None of the nominated loci carries a curated event._")
    else:
        A("| gene | trait | PP4 sQTL | tissues coloc | colocalizing introns | interval match | promotes | tier |")
        A("|---|---|---|---|---|---|---|---|")
        for r in sub.itertuples(index=False):
            A(f"| {r.symbol} | {r.trait} | {r.PP4_sQTL:.3f} | "
              f"{r.n_tissue_sQTL_coloc}/{r.n_tissue} | {r.colocalizing_introns or '—'} | "
              f"{'yes' if r.curated_interval_matched else 'no'} | "
              f"{'**yes**' if r.curated_matched else 'no'} | `{r.tier}` |")
    A("")
    A("## Top nominations by posterior")
    A("")
    A("| gene | trait | PP4 sQTL | PP4 eQTL | introns named | driver tx | tier |")
    A("|---|---|---|---|---|---|---|")
    for r in audit.sort_values("PP4_sQTL", ascending=False).head(30).itertuples(index=False):
        A(f"| {r.symbol} | {r.trait} | {r.PP4_sQTL:.3f} | {r.PP4_eQTL:.3f} | "
          f"{r.n_introns_named} | {r.n_driver_transcripts} | `{r.tier}` |")
    A("")
    A("## References for the curated events")
    A("")
    A("Every interval above is traceable to a publication and to the local file it was "
      "read from. Nothing in the registry is written from recall.")
    A("")
    for r in curated.itertuples(index=False):
        A(f"**{r.gene}** ({r.trait}) — `{r.status}`")
        for src in (r.sources or []):
            pmid = src.get("pmid", "?")
            doi = src.get("doi", "")
            A(f"- PMID [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)"
              + (f" · DOI [{doi}](https://doi.org/{doi})" if doi else ""))
            if src.get("note"):
                A(f"  - {' '.join(str(src['note']).split())}")
            if src.get("supplementary"):
                A(f"  - Supplementary: {' '.join(str(src['supplementary']).split())}")
        for field, label in (("which_event_and_why", "Which event, and why"),
                             ("verification", "Coordinate verification"),
                             ("gtex_phenotype_note", "GTEx testability"),
                             ("scope", "Scope"),
                             ("caution", "Caution")):
            v = getattr(r, field, None)
            if isinstance(v, str) and v.strip():
                A(f"- *{label}:* {' '.join(v.split())}")
        A("")
    A("Article metadata retrieved from PubMed.")
    A("")
    (dest / "LOCUS_EVENT_AUDIT.md").write_text("\n".join(L) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=("audit",), default="audit")
    ap.add_argument("--nominations", choices=("abf", "susie"), default="abf",
                    help="which colocalization layer supplies the nominated loci")
    ap.add_argument("--pp4-call", type=float, default=PP4_CALL)
    ap.add_argument("--sqtl-arm", choices=("representative", "all"),
                    default="representative",
                    help="which sQTL phenotype arm produced the nominations. Only `all` "
                         "puts a non-representative curated event on trial, so only it "
                         "can assign context_distinct_splice_colocalization.")
    ap.add_argument("--max-snps", type=int, default=MAX_SNPS_PRIMARY,
                    help="GWAS SuSiE SNP guard of the coloc.susie results to audit. Anything "
                         "but the primary reads, and writes, a sensitivity arm.")
    args = ap.parse_args(argv)
    run_audit(nominations=args.nominations, call=args.pp4_call,
              sqtl_arm=args.sqtl_arm, max_snps=args.max_snps)


if __name__ == "__main__":
    main()
