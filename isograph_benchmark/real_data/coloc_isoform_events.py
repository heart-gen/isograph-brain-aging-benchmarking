"""Join each colocalized gene's GWAS→sQTL direction to IsoGraph's isoform switch.

`coloc_direction.py` resolves, per colocalized switch gene, the GTEx LeafCutter junction
and the risk-allele-signed direction. This module ties that back to the IsoGraph switch
itself, answering "what is the actual isoform event?" per prioritized gene by joining:

  * GTEx side  — the colocalizing junction (sQTL) + risk-allele direction (coloc_direction).
  * IsoGraph side (interpret_modules, GENCODE v47), taken from the SAME GTEx brain tissue
    the hit colocalized in (tissue-matched, most defensible):
      - structure_switch_pairs : the named switching transcript pair(s) per gene.
      - transcript_polarity_table : per transcript, switch DIRECTION (`r`; opposite signs
        mark the two switch poles), `switch_strength`, and the structural consequence
        (first/last/internal exon change, CDS/UTR change, biotype/coding-status change).
  * GENCODE v47 gtf cache — per-transcript exon coords -> intron set, so a GTEx junction
    can be mapped to the transcript(s) that actually splice it out.
  * BrainSeq independent cohort — where a BrainSeq brain region matches the hit's tissue
    (caudate/hippocampus/dlpfc; caudate_sczd for the disease case) and carries the same
    prioritized gene, we check whether the junction transcript is also an IsoGraph
    switch-pair isoform there. That is out-of-GTEx replication of the isoform switch.

The key added value is the **concordance** flag: does the GTEx sQTL junction map onto a
transcript that IsoGraph flagged as a member of the gene's switching pair *in that tissue*?
If yes, the sQTL and the IsoGraph switch concern the same isoforms — the risk allele shifts
usage toward the junction-containing pole, and we can name the structural event it causes.

Output (<analysis>/coloc/): coloc_isoform_events.parquet + COLOC_ISOFORM_EVENTS.md, plus
a combined parquet across analyses.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import cohort_dir, stage_out
from isograph_benchmark.real_data import gwas_traits as gt
from isograph_benchmark.real_data.qtl_anchoring import _GTEX_TISSUE, _bare

_COLOC_ROOT = stage_out("anchoring.coloc")
_GTEX_ROOT = cohort_dir("gtex")
_BRAINSEQ_ROOT = cohort_dir("brainseq")
_GTF_CACHE = stage_out("tmp",
                 "gencode.v47.primary_assembly.annotation.gtf_cache.parquet")

# GTEx tissue name (as recorded on each coloc hit) -> region artifact dir.
_TISSUE_REGION = {v: k for k, v in _GTEX_TISSUE.items()}

# GTEx tissue -> matching independent BrainSeq cohort region(s). BrainSeq covers caudate,
# hippocampus, and DLPFC (~ frontal cortex); caudate_sczd is the SCZ-diagnosis caudate.
_TISSUE_BRAINSEQ = {
    "Brain_Caudate_basal_ganglia": ["caudate", "caudate_sczd"],
    "Brain_Hippocampus": ["hippocampus"],
    "Brain_Frontal_Cortex_BA9": ["dlpfc"],
    "Brain_Cortex": ["dlpfc"],
}

# LeafCutter intron vs GENCODE exon-boundary conventions can differ by a base; allow a
# small tolerance when matching a junction to a transcript's spliced intron.
_JUNC_TOL = 2

_STRUCT_FLAGS = ["first_exon_changed", "last_exon_changed", "internal_exon_difference",
                 "cds_changed", "utr_changed", "biotype_switch", "coding_status_change"]


# --------------------------------------------------------------------------- GTF
def _transcript_introns(
        genes_bare: set[str]) -> dict[str, dict[str, list[tuple[str, int, int]]]]:
    """bare gene_id -> {transcript_id (versioned): sorted spliced introns}.

    An intron between consecutive exons (exon_i.end, exon_{i+1}.start) is [end+1, start-1]
    in 1-based inclusive coords. Keyed by gene so a junction only matches its own gene's
    transcripts (never another gene's transcript at coincidental coordinates).
    """
    gc = pd.read_parquet(_GTF_CACHE, columns=["transcript_id", "gene_id", "chrom",
                                              "feature", "start", "end"])
    gc = gc[gc["feature"] == "exon"].copy()
    gc["gene_bare"] = _bare(gc["gene_id"])
    gc = gc[gc["gene_bare"].isin(genes_bare)]
    by_gene: dict[str, dict[str, list[tuple[str, int, int]]]] = {}
    for (gene, tx), ex in gc.groupby(["gene_bare", "transcript_id"]):
        ex = ex.sort_values("start")
        chrom = ex["chrom"].iloc[0]
        starts = ex["start"].to_numpy()
        ends = ex["end"].to_numpy()
        tx_introns = [(chrom, int(ends[i]) + 1, int(starts[i + 1]) - 1)
                      for i in range(len(ex) - 1)]
        if tx_introns:
            by_gene.setdefault(gene, {})[tx] = tx_introns
    return by_gene


def _match_junction(chrom: str, start: int, end: int,
                    gene_introns: dict[str, list[tuple[str, int, int]]]) -> list[str]:
    """Transcripts of one gene whose spliced intron equals the junction (± tolerance)."""
    hits = []
    for tx, tx_introns in gene_introns.items():
        for (c, s, e) in tx_introns:
            if c == chrom and abs(s - start) <= _JUNC_TOL and abs(e - end) <= _JUNC_TOL:
                hits.append(tx)
                break
    return hits


def _parse_intron_event(event: str) -> tuple[str, int, int] | None:
    """`chr11:6611818-6611958(-)` -> ('chr11', 6611818, 6611958)."""
    if not isinstance(event, str) or ":" not in event or "-" not in event:
        return None
    try:
        chrom, rest = event.split(":", 1)
        span = rest.split("(", 1)[0]
        s, e = span.split("-")
        return chrom, int(s), int(e)
    except (ValueError, IndexError):
        return None


# ----------------------------------------------------------------- IsoGraph evidence
class _RegionEvidence:
    """IsoGraph switch-pair membership + best-polarity rows for one region artifact dir."""

    def __init__(self, root: Path, region: str):
        self.region = region
        self.members: dict[str, set[str]] = {}        # gene -> switch-pair transcript ids
        self.n_pairs: dict[str, int] = {}             # gene -> #switch pairs
        self.pol: dict[tuple[str, str], pd.Series] = {}  # (gene, tx) -> polarity row
        self.struct: dict[str, pd.Series] = {}        # transcript_id -> structural flags
        mi = root / region / "_m" / "isograph_vae" / "module_interpret"
        sp = mi / "structure_switch_pairs.parquet"
        if sp.exists():
            d = pd.read_parquet(sp)
            d["gene_bare"] = _bare(d["gene_id"])
            for g, sub in d.groupby("gene_bare"):
                self.members[g] = set(sub["transcript_id_1"]) | set(sub["transcript_id_2"])
                self.n_pairs[g] = len(sub)
        sa = mi / "structure_annotations.parquet"
        if sa.exists():
            a = pd.read_parquet(sa).drop_duplicates("transcript_id")
            self.struct = {r.transcript_id: r for r in a.itertuples()}
        pol = []
        for tpt in mi.glob("*/transcript_polarity_table.parquet"):
            pol.append(pd.read_parquet(tpt))
        if pol:
            p = pd.concat(pol, ignore_index=True)
            p["gene_bare"] = _bare(p["gene_id"])
            p = p.sort_values("switch_strength", ascending=False).drop_duplicates(
                ["gene_bare", "transcript_id"])
            self.pol = {(r.gene_bare, r.transcript_id): r for r in p.itertuples()}


def _get_evidence(cache: dict, root: Path, region: str) -> _RegionEvidence:
    key = (str(root), region)
    if key not in cache:
        cache[key] = _RegionEvidence(root, region)
    return cache[key]


def _switch_struct_label(ev: "_RegionEvidence", members: set[str]) -> str:
    """Structural nature of a switch = union of True flags across its pair isoforms."""
    flags: set[str] = set()
    for tx in members:
        sr = ev.struct.get(tx)
        if sr is None:
            continue
        d = sr._asdict()
        for f in _STRUCT_FLAGS:
            if d.get(f) is True:
                flags.add(f.replace("_changed", "").replace("_difference", "")
                          .replace("_", " "))
    return ", ".join(sorted(flags)) if flags else "no annotated structural change"


# --------------------------------------------------------------------------- main run
def run(analysis: str, introns: dict, evidence_cache: dict) -> pd.DataFrame:
    coloc_dir = _COLOC_ROOT / analysis / "coloc"
    dpath = coloc_dir / "coloc_direction.parquet"
    if not dpath.exists():
        print(f"{analysis}: no coloc_direction.parquet; run coloc_direction first.")
        return pd.DataFrame()
    direction = pd.read_parquet(dpath)
    if direction.empty:
        return direction
    trait = direction["trait"].iloc[0]
    case = gt.get(trait).case

    rows = []
    for _, d in direction.iterrows():
        g = d["gene"]
        gene_introns = introns.get(g, {})
        region = _TISSUE_REGION.get(d["tissue"])
        ev = _get_evidence(evidence_cache, _GTEX_ROOT, region) if region else None
        members = ev.members.get(g, set()) if ev else set()

        rec = {
            "analysis": analysis, "gene_source": analysis.rsplit("__", 1)[0],
            "trait": trait, "case": case,
            "gene": g, "gene_name": d.get("gene_name"), "kind": d["kind"],
            "tissue": d["tissue"], "iso_region": region or "",
            "best_rsid": d["best_rsid"], "risk_allele": d.get("risk_allele"),
            "risk_qtl_effect": d.get("risk_qtl_effect"),
            "junction": d.get("intron_event"), "clpp": d.get("clpp"),
            "go_invisible": d.get("go_invisible"),
            "switch_pair": " | ".join(sorted(members)) if members else "",
            "n_switch_pairs": ev.n_pairs.get(g, 0) if ev else 0,
        }

        matched, matched_in_pair, matched_r = [], False, []
        parsed = _parse_intron_event(d.get("intron_event")) if d["kind"] == "sQTL" else None
        if parsed is not None:
            chrom, s, e = parsed
            matched = _match_junction(chrom, s, e, gene_introns)
            for tx in matched:
                matched_in_pair = matched_in_pair or (tx in members)
                pr = ev.pol.get((g, tx)) if ev else None
                if pr is not None:
                    matched_r.append(float(pr.r))
        rec["junction_transcripts"] = ",".join(matched)
        rec["junction_in_switch_pair"] = matched_in_pair
        rec["junction_transcript_polarity_r"] = float(np.mean(matched_r)) if matched_r else np.nan
        # Structural nature of the switch = union of changes across its pair isoforms
        # (flags in structure_annotations are per-transcript vs a reference, so the switch
        # is described by pooling the True flags over the switching pair members).
        rec["structural_consequence"] = (
            _switch_struct_label(ev, members) if (ev and matched_in_pair) else "")
        rec["concordant"] = (matched_in_pair if (d["kind"] == "sQTL" and matched)
                             else pd.NA)

        # BrainSeq independent-cohort replication: same gene, matching region, junction
        # transcript is also an IsoGraph switch-pair isoform out of GTEx.
        bs_region, bs_members, bs_rep = "", set(), pd.NA
        for cand in _TISSUE_BRAINSEQ.get(d["tissue"], []):
            # for disease hits prefer the diagnosis cohort; skip caudate_sczd for aging.
            if cand == "caudate_sczd" and case != "disease":
                continue
            if cand != "caudate_sczd" and case == "disease" and cand == "caudate":
                pass  # aging caudate still a valid independent cohort for a disease hit
            bev = _get_evidence(evidence_cache, _BRAINSEQ_ROOT, cand)
            m = bev.members.get(g, set())
            if m:
                bs_region, bs_members = cand, m
                if matched:
                    bs_rep = any(tx in m for tx in matched)
                break
        rec["brainseq_region"] = bs_region
        rec["brainseq_switch_pair"] = " | ".join(sorted(bs_members)) if bs_members else ""
        rec["brainseq_replicates_switch"] = bs_rep

        rec["resolved_event"] = _event_sentence(d, rec)
        rows.append(rec)

    out = pd.DataFrame(rows).sort_values(["kind", "clpp"], ascending=[True, False])
    out.to_parquet(coloc_dir / "coloc_isoform_events.parquet", index=False)
    _write_report(coloc_dir, analysis, out)
    n_conc = int((out["concordant"] == True).sum())
    n_rep = int((out["brainseq_replicates_switch"] == True).sum())
    print(f"{analysis}: {len(out)} events; GTEx junction↔switch concordant: {n_conc}; "
          f"BrainSeq-replicated: {n_rep}")
    return out


def _event_sentence(d: pd.Series, rec: dict) -> str:
    eff = rec["risk_qtl_effect"]
    if pd.isna(eff):
        return "direction unresolved"
    ra, arrow = rec["risk_allele"], ("increases" if eff > 0 else "decreases")
    if d["kind"] == "eQTL":
        return f"risk allele {ra} {arrow} {rec['gene_name']} expression (gene-level; no intron)"
    tx = rec["junction_transcripts"] or "unmapped transcript"
    conc = ("matches an IsoGraph switch-pair isoform" if rec["junction_in_switch_pair"]
            else "not in the IsoGraph switch pair for this tissue")
    rep = ""
    if rec["brainseq_replicates_switch"] is True:
        rep = f"; replicated in BrainSeq {rec['brainseq_region']}"
    struct = rec["structural_consequence"] or "structural change n/a"
    return (f"risk allele {ra} {arrow} usage of junction {rec['junction']} "
            f"({tx}; {conc}{rep}) — {struct}")


def _write_report(coloc_dir: Path, analysis: str, out: pd.DataFrame) -> None:
    sq = out[out["kind"] == "sQTL"]
    conc = sq[sq["concordant"] == True]
    rep = sq[sq["brainseq_replicates_switch"] == True]
    lines = [
        f"# Resolved isoform events — {analysis}", "",
        "Each colocalized switch gene: the GTEx risk-allele direction joined to the "
        "IsoGraph switch (structure_switch_pairs + transcript_polarity) from the SAME GTEx "
        "tissue, GENCODE v47. `concordant` = the sQTL junction maps onto an IsoGraph "
        "switch-pair isoform in that tissue; `BrainSeq` = the same junction transcript is "
        "an IsoGraph switch-pair isoform in the matching independent BrainSeq cohort.", "",
        f"- total colocalized events: **{len(out)}**",
        f"- sQTL events with a mapped junction transcript: "
        f"**{int((sq['junction_transcripts'] != '').sum())}/{len(sq)}**",
        f"- sQTL junction ↔ IsoGraph switch-pair concordant (GTEx): **{len(conc)}**",
        f"- also replicated in an independent BrainSeq cohort: **{len(rep)}**", "",
        "| gene | trait | kind | tissue | junction | risk→ | GTEx | BrainSeq | structural event |",
        "|------|-------|------|--------|----------|-------|------|----------|------------------|",
    ]
    for _, r in out.head(30).iterrows():
        sign = ("↑" if (pd.notna(r["risk_qtl_effect"]) and r["risk_qtl_effect"] > 0)
                else "↓" if pd.notna(r["risk_qtl_effect"]) else "?")
        c = "yes" if r["concordant"] is True else ("no" if r["concordant"] is False else "—")
        b = ("yes" if r["brainseq_replicates_switch"] is True
             else "no" if r["brainseq_replicates_switch"] is False else "—")
        lines.append(
            f"| {r['gene_name']} | {r['trait'].upper()} | {r['kind']} | {r['iso_region']} | "
            f"{r['junction'] or '(gene-level)'} | {sign} | {c} | {b} | "
            f"{r['structural_consequence'] or '—'} |")
    (coloc_dir / "COLOC_ISOFORM_EVENTS.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Resolve the isoform event per colocalized gene (join sQTL direction "
                    "to the tissue-matched IsoGraph switch + BrainSeq replication).")
    p.add_argument("--analysis", nargs="*", default=None,
                   help="analysis dir(s); default: all with coloc_direction.parquet")
    args = p.parse_args()
    analyses = args.analysis or sorted(
        d.name for d in _COLOC_ROOT.iterdir()
        if (d / "coloc" / "coloc_direction.parquet").exists())

    all_genes: set[str] = set()
    for a in analyses:
        dp = _COLOC_ROOT / a / "coloc" / "coloc_direction.parquet"
        all_genes |= set(pd.read_parquet(dp, columns=["gene"])["gene"])
    introns = _transcript_introns(all_genes)
    evidence_cache: dict = {}

    combined = []
    for a in analyses:
        df = run(a, introns, evidence_cache)
        if not df.empty:
            combined.append(df)
    if combined:
        allout = pd.concat(combined, ignore_index=True)
        allout.to_parquet(_COLOC_ROOT / "coloc_isoform_events_combined.parquet", index=False)
        print(f"combined -> {_COLOC_ROOT/'coloc_isoform_events_combined.parquet'} "
              f"({len(allout)} events)")
        _write_meta(allout)


def _write_meta(allout: pd.DataFrame) -> None:
    """Cross-trait rollup: per trait, colocalized events and how many resolve to an
    IsoGraph switch-pair isoform (GTEx-concordant) or replicate in a BrainSeq cohort."""
    rows = []
    # key on (gene_source, trait): the aging switch layer projected onto SCZ (aging__scz)
    # and the disease switch layer on SCZ (brainseq-sczd__scz) share trait=scz but are
    # distinct anchoring analyses and must stay separate rows.
    for (source, trait), sub in allout.groupby(["gene_source", "trait"]):
        sq = sub[sub["kind"] == "sQTL"]
        conc = sq[sq["concordant"] == True]
        rows.append({
            "gene_source": source,
            "trait": trait.upper(),
            "case": sub["case"].iloc[0],
            "events": len(sub),
            "sQTL_events": len(sq),
            "junction_mapped": int((sq["junction_transcripts"] != "").sum()),
            "GTEx_concordant": len(conc),
            "concordant_go_invisible": int((conc["go_invisible"] == True).sum()),
            "brainseq_replicated": int((sub["brainseq_replicates_switch"] == True).sum()),
        })
    meta = pd.DataFrame(rows).sort_values(["case", "gene_source", "trait"])
    meta.to_parquet(_COLOC_ROOT / "coloc_isoform_events_meta.parquet", index=False)
    lines = [
        "# Resolved isoform events — cross-trait rollup", "",
        "Per anchoring analysis (switch layer × trait): colocalized events, how many sQTL "
        "junctions map to a transcript, how many of those map onto an IsoGraph switch-pair "
        "isoform in the same GTEx tissue (`GTEx_concordant`; `_go_invisible` = in a "
        "GO-invisible switch module), and how many replicate the switch in an independent "
        "BrainSeq cohort.", "",
        "| switch layer | trait | case | events | sQTL | junction mapped | GTEx concordant | GO-invisible | BrainSeq replicated |",
        "|--------------|-------|------|--------|------|-----------------|-----------------|--------------|---------------------|",
    ]
    for _, r in meta.iterrows():
        lines.append(
            f"| {r['gene_source']} | {r['trait']} | {r['case']} | {r['events']} | "
            f"{r['sQTL_events']} | {r['junction_mapped']} | {r['GTEx_concordant']} | "
            f"{r['concordant_go_invisible']} | {r['brainseq_replicated']} |")
    (_COLOC_ROOT / "COLOC_ISOFORM_EVENTS_META.md").write_text("\n".join(lines) + "\n")
    print(f"meta -> {_COLOC_ROOT/'COLOC_ISOFORM_EVENTS_META.md'}")


if __name__ == "__main__":
    main()
