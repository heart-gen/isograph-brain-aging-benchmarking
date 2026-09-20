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

`--layer signal` resolves the signal-level (`coloc.susie`, all-introns arm) nominations
instead, through the same resolver, and writes beside them under
`coloc_signal_susie/all_introns/`. It is a second event table, not a replacement: the CLPP
table stays what the committed displays read, and the signal table carries no direction.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from isograph_benchmark.paths import cohort_dir, stage_out
from isograph_benchmark.real_data import gwas_traits as gt
from isograph_benchmark.real_data.coloc_signal_susie import (
    prior_robustness,
    results_dir,
    signal_root,
)
from isograph_benchmark.real_data.locus_event_audit import (
    CONTEXT_DISTINCT_ARM,
    PP4_CALL,
    _parse_gtex_intron,
    load_nominations,
)
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


def _flag_true(v) -> bool:
    """True only for a definite True, whatever container the flag arrived in.

    `structure_annotations` stores these columns as the pandas nullable `boolean`
    dtype, so `itertuples` hands back `numpy.bool_`. `numpy.bool_(True) is True` is
    False, so an identity test silently dropped every flag and labelled all 76
    concordant events "no annotated structural change". A missing value must still be
    rejected without being coerced, because `bool(pd.NA)` raises.
    """
    return bool(pd.notna(v) and bool(v))


def _switch_struct_label(ev: "_RegionEvidence", members: set[str]) -> str:
    """Structural nature of a switch = union of True flags across its pair isoforms."""
    flags: set[str] = set()
    for tx in members:
        sr = ev.struct.get(tx)
        if sr is None:
            continue
        d = sr._asdict()
        for f in _STRUCT_FLAGS:
            if _flag_true(d.get(f)):
                flags.add(f.replace("_changed", "").replace("_difference", "")
                          .replace("_", " "))
    return ", ".join(sorted(flags)) if flags else "no annotated structural change"


# --------------------------------------------------------------------------- main run
def _resolve_event(d, analysis: str, trait: str, case: str, introns: dict,
                   evidence_cache: dict) -> dict:
    """One colocalized (gene, tissue, junction) -> the IsoGraph switch it lands on.

    Both coloc layers resolve through here, so a CLPP event and a signal-level event are
    matched to transcripts, to tissue-matched switch evidence and to BrainSeq by identical
    code. `d` carries gene, gene_name, kind, tissue and intron_event, plus whatever
    direction and CLPP fields its layer has; an absent field resolves to missing.
    """
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
        "best_rsid": d.get("best_rsid"), "risk_allele": d.get("risk_allele"),
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
    return rec


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
        rec = _resolve_event(d, analysis, trait, case, introns, evidence_cache)
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


# ----------------------------------------------------------------- signal-level layer
def signal_events_path() -> Path:
    """The signal-level event table: beside the all-introns nominations it resolves, and
    never on top of the CLPP layer's combined table, which the committed displays read."""
    return results_dir(signal_root("switch"), "all") / "coloc_isoform_events.parquet"


def signal_event_cells(nom: pd.DataFrame, cells: pd.DataFrame, srep: pd.DataFrame,
                       call: float = PP4_CALL) -> pd.DataFrame:
    """The per-tissue sQTL colocalizations behind each signal-level nomination.

    This is the selection the locus event audit tiers on, so both name the same introns:
    a nominated (analysis, trait, locus, gene), each tissue where its sQTL PP4 reaches the
    call, and the intron that tissue's call used. A `coloc.abf` fallback cell has no fitted
    intron, so it names GTEx's representative one and keeps `estimator_sQTL = abf` -- it
    resolves an event, but it tested only that intron.
    """
    keys = ["analysis", "trait", "LOCUS_ID", "gene"]
    c = cells[cells["PP4_sQTL"] >= call].merge(nom[keys].drop_duplicates(), on=keys,
                                              how="inner")
    rep = (srep[["gene", "tissue", "phenotype_id"]].drop_duplicates(["gene", "tissue"])
           .rename(columns={"phenotype_id": "_rep"}))
    c = c.merge(rep, on=["gene", "tissue"], how="left")
    pid = c["phenotype_id"] if "phenotype_id" in c.columns else pd.Series(None, index=c.index)
    c["phenotype_id"] = pid.where(pid.notna(), c["_rep"])
    c = c.drop(columns="_rep")
    parsed = c["phenotype_id"].map(_parse_gtex_intron)
    c = c[parsed.notna()].copy()

    def _event(pid: str) -> str:
        chrom, s, e = _parse_gtex_intron(pid)
        f = pid.split(":")
        strand = f[3].rsplit("_", 1)[-1] if len(f) > 3 else ""
        return f"{chrom}:{s}-{e}({strand})" if strand in ("+", "-") else f"{chrom}:{s}-{e}"

    c["intron_event"] = c["phenotype_id"].map(_event)
    return c.reset_index(drop=True)


def _is_true(s: pd.Series) -> pd.Series:
    return s.map(lambda v: v is True)


def _signal_event_sentence(rec: dict) -> str:
    tx = rec["junction_transcripts"] or "unmapped transcript"
    conc = ("matches an IsoGraph switch-pair isoform" if rec["junction_in_switch_pair"]
            else "not in the IsoGraph switch pair for this tissue")
    rep = (f"; replicated in BrainSeq {rec['brainseq_region']}"
           if rec["brainseq_replicates_switch"] is True else "")
    struct = rec["structural_consequence"] or "structural change n/a"
    return (f"junction {rec['junction']} colocalizes (PP4 {rec['PP4_sQTL']:.2f}, "
            f"{rec['estimator']}; {tx}; {conc}{rep}) — {struct}; direction not resolved")


# Genes reviewed in the long-read arm even when the colocalizing tissue carries no IsoGraph
# switch pair for them (user decision, 2026-09-11). UNC13A's ALS signal colocalizes only in
# the two cerebellar tissues, where IsoGraph called no UNC13A switch, yet the transcripts that
# splice its colocalizing intron are switch-pair members in other GTEx brain regions,
# including frontal cortex BA9 -- the long-read tissue. The exception changes what the test
# means, not whether it runs: the pair is then a switch defined in another region, so these
# rows are flagged `cross_tissue_exception`, never counted as tissue-concordant, and scored
# apart from the set-level statistic.
CROSS_TISSUE_EXCEPTIONS: tuple[str, ...] = ("UNC13A",)


def cross_tissue_switch_pairs(junction_tx: set[str], gene: str,
                              members_by_region: dict[str, dict[str, set[str]]]
                              ) -> tuple[list[str], set[str]]:
    """Regions whose switch pair for `gene` contains a junction transcript, and the union of
    those pairs' members. Pure: `members_by_region` maps region -> gene -> pair members."""
    regions, members = [], set()
    for region in sorted(members_by_region):
        m = members_by_region[region].get(gene, set())
        if m & junction_tx:
            regions.append(region)
            members |= m
    return regions, members


def _cross_tissue_fields(rec: dict, symbol, exceptions: set[str], cache: dict) -> dict:
    out = {"switch_pair_scope": "tissue_matched" if rec["junction_in_switch_pair"] else "none",
           "cross_tissue_regions": "", "cross_tissue_switch_pair": "",
           "cross_tissue_in_switch_pair": False}
    if (rec["junction_in_switch_pair"] or str(symbol) not in exceptions
            or not rec["junction_transcripts"]):
        return out
    by_region = {r: _get_evidence(cache, _GTEX_ROOT, r).members
                 for r in sorted(set(_TISSUE_REGION.values()))}
    regions, members = cross_tissue_switch_pairs(
        set(rec["junction_transcripts"].split(",")), rec["gene"], by_region)
    if regions:
        out.update(switch_pair_scope="cross_tissue_exception",
                   cross_tissue_regions=";".join(regions),
                   cross_tissue_switch_pair=" | ".join(sorted(members)),
                   cross_tissue_in_switch_pair=True)
    return out


def run_signal(call: float = PP4_CALL,
               cross_tissue_exceptions: tuple[str, ...] = CROSS_TISSUE_EXCEPTIONS
               ) -> pd.DataFrame:
    nom, cells = load_nominations("susie", call=call, sqtl_arm=CONTEXT_DISTINCT_ARM)
    srep_f = stage_out("anchoring.coloc_modality") / "sqtl_representative.parquet"
    srep = (pd.read_parquet(srep_f) if srep_f.exists()
            else pd.DataFrame(columns=["gene", "tissue", "phenotype_id"]))
    ec = signal_event_cells(nom, cells, srep, call)
    if ec.empty:
        raise SystemExit("no signal-level sQTL colocalizations to resolve")
    introns = _transcript_introns(set(ec["gene"]))
    cache: dict = {}
    exceptions = {str(g) for g in cross_tissue_exceptions}

    rows = []
    for d in ec.to_dict("records"):
        trait = str(d["trait"])
        rec = _resolve_event({**d, "gene_name": d.get("symbol"), "kind": "sQTL"},
                             str(d["analysis"]), trait, gt.get(trait).case, introns, cache)
        # A CLPP column here would be a posterior this layer never computed, and no signed
        # direction exists at signal level (`coloc_direction` is the CLPP layer's).
        del rec["clpp"]
        rec.update(_cross_tissue_fields(rec, d.get("symbol"), exceptions, cache))
        p12 = d.get("p12_min_call")
        p12 = float(p12) if pd.notna(p12) else np.nan
        rec.update({
            "LOCUS_ID": d["LOCUS_ID"], "phenotype_id": d["phenotype_id"],
            "PP4_sQTL": float(d["PP4_sQTL"]),
            "PP4_eQTL": float(d["PP4_eQTL"]) if pd.notna(d.get("PP4_eQTL")) else np.nan,
            "estimator": d.get("estimator_sQTL"),
            "p12_min_call": p12,
            "prior_robustness": prior_robustness(p12),
            "fallback_reason": d.get("fallback_reason"),
            "coloc_layer": "signal",
        })
        rec["resolved_event"] = _signal_event_sentence(rec)
        rows.append(rec)

    out = (pd.DataFrame(rows)
           .sort_values(["trait", "gene_name", "PP4_sQTL"], ascending=[True, True, False])
           .reset_index(drop=True))
    dest = signal_events_path()
    out.to_parquet(dest, index=False)
    _write_signal_report(dest.parent, out, call)
    conc = out[_is_true(out["junction_in_switch_pair"])]
    print(f"signal layer: {len(out)} events over "
          f"{len(out[['trait', 'gene']].drop_duplicates())} gene x trait nominations; "
          f"on a switch-pair isoform: {len(conc)} events, "
          f"{conc['gene'].nunique()} genes -> {dest}")
    return out


def _write_signal_report(dest: Path, out: pd.DataFrame, call: float) -> None:
    conc = out[_is_true(out["junction_in_switch_pair"])]
    n_nom = len(out[["trait", "gene"]].drop_duplicates())
    lines = [
        "# Resolved isoform events — signal-level colocalization (all-introns arm)", "",
        f"Every tissue where a `coloc.susie` nomination's sQTL reaches PP4 >= {call}, the "
        "intron that call used, and the IsoGraph switch it lands on in the same GTEx tissue. "
        "Junction matching, switch evidence and the BrainSeq check are the same code as the "
        "CLPP layer (`coloc/coloc_isoform_events_combined.parquet`), which this table sits "
        "beside and does not replace.", "",
        "Three differences from the CLPP layer, each limiting what a row can say:", "",
        "- **No direction.** The risk-allele-signed effect is resolved only for the CLPP "
        "layer (`coloc_direction`); `risk_qtl_effect` is empty here rather than borrowed.",
        "- **`estimator = abf` rows name GTEx's representative intron**, the only one the "
        "fallback scored. They resolve an event but never tested the gene's other introns.",
        "- **One intron per tissue:** the intron behind that tissue's call. Other introns "
        "of the gene that also colocalize are in `signal_pairs.parquet`, not here.", "",
        f"- events (nomination x tissue): **{len(out)}** over **{n_nom}** gene x trait "
        "nominations",
        f"- signal-level (`susie`): **{int((out['estimator'] == 'susie').sum())}**; "
        f"abf fallback: **{int((out['estimator'] == 'abf').sum())}**",
        f"- junction mapped to a GENCODE v47 transcript: "
        f"**{int((out['junction_transcripts'] != '').sum())}/{len(out)}**",
        f"- junction on an IsoGraph switch-pair isoform in the same tissue: "
        f"**{len(conc)}** events, **{len(conc[['trait', 'gene']].drop_duplicates())}** "
        "gene x trait",
        f"- also a switch-pair isoform in the matching BrainSeq region: "
        f"**{int(_is_true(out['brainseq_replicates_switch']).sum())}**",
        f"- cross-tissue exceptions (switch pair taken from another region; flagged, never "
        f"counted as tissue-concordant, scored apart in the long-read arm): "
        f"**{int(_is_true(out['cross_tissue_in_switch_pair']).sum())}** events "
        f"({', '.join(sorted(out.loc[_is_true(out['cross_tissue_in_switch_pair']), 'gene_name'].astype(str).unique())) or 'none'})",
        "",
        "| gene | trait | tissue | junction | PP4 | estimator | prior | transcripts | "
        "switch pair | BrainSeq | structural event |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in out.itertuples(index=False):
        c = ("yes" if r.concordant is True else "no" if r.concordant is False else "—")
        if r.cross_tissue_in_switch_pair is True:
            c = f"cross-tissue exception ({r.cross_tissue_regions})"
        b = ("yes" if r.brainseq_replicates_switch is True
             else "no" if r.brainseq_replicates_switch is False else "—")
        lines.append(
            f"| {r.gene_name} | {r.trait.upper()} | {r.iso_region} | {r.junction} | "
            f"{r.PP4_sQTL:.3f} | {r.estimator} | {r.prior_robustness} | "
            f"{r.junction_transcripts or '—'} | {c} | {b} | {r.structural_consequence or '—'} |")
    (dest / "COLOC_ISOFORM_EVENTS.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    p = argparse.ArgumentParser(
        description="Resolve the isoform event per colocalized gene (join sQTL direction "
                    "to the tissue-matched IsoGraph switch + BrainSeq replication).")
    p.add_argument("--layer", choices=("clpp", "signal"), default="clpp",
                   help="clpp: eCAVIAR coloc_direction events (default; the layer the "
                        "committed displays read). signal: coloc.susie all-introns "
                        "nominations, written beside them")
    p.add_argument("--analysis", nargs="*", default=None,
                   help="analysis dir(s); default: all with coloc_direction.parquet")
    p.add_argument("--cross-tissue-exception", nargs="*", metavar="SYMBOL",
                   default=list(CROSS_TISSUE_EXCEPTIONS),
                   help="signal layer: genes whose switch pair may come from a region other "
                        "than the colocalizing tissue (flagged, scored apart; default: "
                        f"{' '.join(CROSS_TISSUE_EXCEPTIONS)}). Pass with no symbols to disable")
    args = p.parse_args()
    if args.layer == "signal":
        if args.analysis:
            raise SystemExit("--analysis selects CLPP analysis dirs; the signal layer "
                             "resolves every signal-level nomination")
        run_signal(cross_tissue_exceptions=tuple(args.cross_tissue_exception))
        return
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
