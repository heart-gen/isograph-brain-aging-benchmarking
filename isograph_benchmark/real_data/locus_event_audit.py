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
    disease_locus_splice_linked coloc holds and an event is named, but not that event
    novel_splice_led_candidate  coloc holds, event named, no curated event at the locus
    not_resolved                coloc does not hold, or no intron could be named

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

_EVENTS_YAML = rel("configs", "known_splice_events.yaml")
_GTF_CACHE = stage_out("tmp", "gencode.v47.primary_assembly.annotation.gtf_cache.parquet")

# Pre-specified before the audit was run.
PP4_CALL = 0.80
TIERS = ("known_mechanism_recovered", "disease_locus_splice_linked",
         "novel_splice_led_candidate", "not_resolved")
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


def out_dir(nominations: str) -> Path:
    return ensure_dir(stage_out("anchoring", "locus_event_audit", nominations))


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
# Nominations
# --------------------------------------------------------------------------- #
def load_nominations(source: str, call: float = PP4_CALL) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(gene-level nominations, per-tissue cells) for a nomination source."""
    if source == "abf":
        base = stage_out("anchoring.coloc_modality")
    elif source == "susie":
        base = stage_out("anchoring.coloc_signal")
    else:
        raise SystemExit(f"unknown nomination source {source!r}")
    gf, cf = base / "genes.parquet", base / "cells.parquet"
    for f in (gf, cf):
        if not f.exists():
            raise SystemExit(f"missing {f}; run that layer's --stage meta first")
    genes = pd.read_parquet(gf)
    cells = pd.read_parquet(cf)
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


def assign_tier(introns: list, cur: pd.DataFrame) -> dict:
    """Tier a locus from its colocalizing introns and the curated events at that locus.

    Two INDEPENDENT gates guard the top tier and a locus must clear both:
      * `status` must be `reviewed`  -- the coordinates were verified in this repo.
      * `evidence_class` must be functional or coloc-based -- a TWAS association does not
        establish a shared causal variant, so it may never certify a mechanism.
    `curated_interval_matched` is reported separately so that a coordinate match capped by
    the evidence gate stays visible instead of looking like no match at all.
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

    if not introns:
        tier = "not_resolved"
    elif matched:
        tier = "known_mechanism_recovered"
    elif len(cur):
        tier = "disease_locus_splice_linked"
    else:
        tier = "novel_splice_led_candidate"
    return {"tier": tier, "curated_matched": matched,
            "curated_interval_matched": interval_matched,
            "curated_status": matched_status,
            "curated_evidence_class": matched_evidence}


def run_audit(nominations: str = "abf", call: float = PP4_CALL) -> Path:
    dest = out_dir(nominations)
    curated = load_curated_events()
    nom, cells = load_nominations(nominations, call=call)
    gtf = _load_gtf()

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
        c = cells[(cells["gene"] == gene) & (cells["trait"] == trait)
                  & (cells["PP4_sQTL"] >= call)]
        introns, tissues, tx_all = [], [], set()
        for cr in c.itertuples(index=False):
            tis = str(cr.tissue)
            pid = (str(cr.phenotype_id) if "phenotype_id" in c.columns
                   and pd.notna(getattr(cr, "phenotype_id", None))
                   else _sqtl_phenotype_for(gene, tis, srep))
            parsed = _parse_gtex_intron(pid) if pid else None
            if parsed is None:
                continue
            chrom, istart, iend = parsed
            tissues.append(tis)
            introns.append((chrom, istart, iend, pid))
            tx_all.update(transcripts_spliced_by(chrom, istart, iend, gene, gtf))

        # Link 4 and the curated match. Extracted into `assign_tier` so the rule that
        # decides whether a locus may be called a recovered mechanism is testable on its
        # own, rather than only reachable by running the whole audit.
        cur = curated[(curated["gene"] == sym) & (curated["trait"] == trait)]
        verdict = assign_tier(introns, cur)
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
            "tier": tier,
        })

    audit = pd.DataFrame(rows).sort_values(["tier", "PP4_sQTL"], ascending=[True, False])
    audit.to_parquet(dest / "locus_event_audit.parquet", index=False)
    _write_report(dest, audit, curated, nominations, call)
    print(f"\n  tiers: {audit['tier'].value_counts().to_dict()}")
    print(f"  wrote {dest}")
    return dest


def _write_report(dest: Path, audit: pd.DataFrame, curated: pd.DataFrame,
                  nominations: str, call: float) -> None:
    L: list[str] = []
    A = L.append
    A("# Locus event audit")
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
    A("## Tier counts")
    A("")
    for t in TIERS:
        A(f"- `{t}`: {int((audit['tier'] == t).sum())}")
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
    args = ap.parse_args(argv)
    run_audit(nominations=args.nominations, call=args.pp4_call)


if __name__ == "__main__":
    main()
