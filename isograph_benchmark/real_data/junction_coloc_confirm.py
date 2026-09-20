"""Short-read junction confirmation of the genetically anchored switch pairs.

Every gene whose colocalizing sQTL junction lands in the tissue-matched IsoGraph switch
pair is tested here -- 30 genes over 76 events on the 2026-09-19 re-run, not the two the
arm was originally pointed at. Until 2026-09-20 ``--genes`` defaulted to ``SNCA CTSH``,
so when those two stopped being the display item the arm reported "No testable targets"
and said nothing about the other 28 concordant genes.

The long-read check (``switch_orthogonal_confirm``) could not confirm the anchored
isoforms: SNCA's sits at 0.29% of the gene's ONT output in n = 12 samples, a failure
confounded with the assay, which measures whole transcripts. This CLI runs the test the
long-read data could not: the same junctions in BrainSEQ short-read, which measures **the
junction the sQTL actually tags**, at roughly 40x the sample size (hippocampus n = 452,
DLPFC n = 500, caudate n = 487).

What the panel can and cannot reach
-----------------------------------
BrainSEQ sequences caudate, DLPFC and hippocampus. GTEx brain covers thirteen tissues, and
**58 of the 76 concordant events are cerebellum or cerebellar hemisphere**, which BrainSEQ
does not sequence at all. Those targets are reported with verdict
``no_matched_brainseq_region`` rather than dropped: a locus this panel cannot measure must
not look like a locus that was never nominated. The ceiling on how much of the concordant
set this arm can ever confirm is therefore set by GTEx's tissue-wise power -- cerebellum has
the largest brain sample size -- not by the genes.

This is the NARROWEST of the three junction-level resources in the repository, and a
``junction_not_measured`` or ``no_matched_brainseq_region`` verdict here is not the last
word on a gene. Over the 30 concordant genes (2026-09-20):

  this arm (LIBD PSI events)            9 genes reach a region, 6 have the junction measured
  ONT long-read (longread_switch_confirm)          15 genes ``confirmed`` at gene level
  BrainSEQ junction recount (ase_junction_switch)  18 genes with >=1 testable switch pair
  none of the three                                 8 genes

PRDM2 is the case in point: this arm returns ``junction_not_measured`` because the LIBD
event catalogue stops at ~13,787,067, while ONT detects both transcripts of the drawn pair
(usage Spearman -0.448, ``switch_like``) and the allele-aware recount measures that exact
junction over 292 donors. Report the three together; do not cite this arm's silence alone.

What is tested
--------------
A colocalizing sQTL says the risk allele moves usage of one junction. For that to be a
*switch* -- and for the locus to carry a main figure -- the junction's isoform must
actually be used, and its usage must trade off against the alternative form. Two tests,
picked by what the PSI table can measure for each locus:

``direct`` (primary, used whenever available)
    A single LIBD PSI event whose two arms contrast the anchored junction against an
    alternative form. SNCA's AF event
    ``AF:chr4:89835692-89836127:89836213:89835692-89838252:89838315:-`` is exactly the
    Fig 4A contrast -- anchored proximal first exon versus canonical distal first exon --
    and CTSH's ``A3:chr15:78937496-78939140:78937423-78939140:-`` contrasts the anchored
    acceptor against the alternative one. For such an event PSI *is* the switch ratio, so
    there is no compositional-closure confound at all: the two arms are the whole of the
    event by construction.

    The statistic is deliberately convention-free. PSI orientation (which arm is the
    numerator) is not documented for these tables, so instead of assuming it we report
    **minor-form usage**, ``min(median PSI, 1 - median PSI)``. A switch is usable only if
    *both* forms carry non-trivial usage; if PSI sits at 0.99 or 0.01 then one form is
    effectively absent whichever way round it is. This is the direct short-read answer to
    the long-read failure, where SNCA's anchored isoform was 0.29% of gene output.

``paired`` (fallback)
    Where the anchored and competing junctions are measured by *different* events, PSI
    values are compared across individuals and the pair must **anti-correlate**. Here the
    compositional-closure baseline matters: two PSI events of one gene are pushed towards
    negative correlation for arithmetic reasons alone, so rho is tested against the
    empirical distribution of within-gene event-pair correlations (``within_gene``) and of
    pairs from genes matched on event count (``matched_gene``) -- never against zero. Both
    PSI vectors are residualized on the standard aging covariates first.

Decision rule (pre-registered, from ORTHOGONAL_CONFIRMATION.md)
--------------------------------------------------------------
If the junctions validate -- anti-correlated beyond the closure baseline at the stated
alpha -- the short-read result is reported as the orthogonal confirmation and the
long-read failure is cited as an assay limitation. If they do not, the locus stays off any
main figure and the set-level result stands alone. The rule is applied by this CLI and
written to disk as ``verdict``; it is not re-decided after seeing the numbers.

The rule was written for a two-gene run, so applied across the whole panel its alpha is an
uncorrected family. ``verdict`` is therefore kept exactly as pre-registered and the
BH-corrected call is added beside it as ``verdict_bh`` (``q_within_gene``,
``q_matched_gene``), which is what the panel-wide count in the report uses. ``direct`` rows
are not in that family: their verdict is a usage threshold, not a test.

The display-item decision rule in ``_REFERENCE_ARM`` stays gene-specific by design -- it
asks whether the exact contrast a figure panel draws holds up, which only means anything
for a gene that has a panel. Genes without one are reported, not decided.

Known limitation of the usage threshold (PI, 2026-09-20)
--------------------------------------------------------
``min_usage`` treats a rare minor form as a failure. BrainSEQ and GTEx are neurotypical
tissue, so an isoform that matters in disease is often a MINORITY form here precisely
because the disease is what raises it -- and this threshold will call exactly those
``not_validated``. It was pre-registered against a different failure (SNCA's anchored
isoform at 0.29% of ONT gene output, where the form was effectively absent rather than
merely rare), so it is left in place as the pre-registered rule and NOT silently rewritten;
changing a pre-registered criterion after seeing the numbers is the PI's call.

Read a ``not_validated`` here as "this form is rare in normal brain", not as "this switch
is not real". ``08_integration/_m/anchored_gene_summary/`` reports the same junctions with
detection counts (how many donors carry BOTH forms) instead of an abundance gate, which is
the criterion that survives the argument above.

Outputs land in ``06_switch_mechanism/_m/junction_coloc_confirm/``.
"""
from __future__ import annotations

import argparse
import json
import re

import numpy as np
import pandas as pd
from scipy import stats

from isograph_benchmark.paths import ensure_dir, stage_out
from isograph_benchmark.real_data.validate_switch_splicing import (
    _bundle_path,
    _covariates,
    _load_splice,
    _parse_junction,
)

# Coordinate tolerance when matching a coloc junction to a PSI event, in bp. Kept
# identical to validate_switch_splicing._COORD_MATCH_TOL so the two CLIs agree on what
# "the same junction" means.
_COORD_TOL = 5

# Every "start-end" pair in an event_info string.
#
# NB this deliberately does NOT reuse validate_switch_splicing's
# ``_JUNC_RE = (chr[\w]+):(\d+)-(\d+)``. A LIBD event_info names the chromosome ONCE and
# then lists both arms, e.g.
#     A3:chr15:78937496-78939140:78937423-78939140:-
# so a chr-anchored pattern matches only the FIRST arm and the competing arm of every
# event is invisible. That is why CTSH's anchored junction chr15:78937423-78939140 --
# which is the second arm of that A3 event -- returned "junction_not_measured" until this
# parser replaced it. The chromosome is checked separately against the event's `chrom`
# column, so dropping it from the pattern loses no specificity.
_PAIR_RE = re.compile(r"(\d+)-(\d+)")


def _event_coords(event_info: str) -> list[int]:
    """All coordinates named by an event_info string, both arms included."""
    return [int(v) for pair in _PAIR_RE.findall(str(event_info)) for v in pair]


def _is_reference(gene: str, event_info: str) -> bool:
    """Does this event contrast the anchored junction against the arm the figure claims?"""
    ref = _REFERENCE_ARM.get(str(gene))
    if ref is None:
        return False
    return any(abs(c - ref) <= _COORD_TOL for c in _event_coords(event_info))

# The GTEx tissue a colocalization was found in -> the BrainSEQ region(s) that measure the
# same anatomy, with how close the match is. BrainSEQ has three regions only -- caudate,
# DLPFC and hippocampus -- so most GTEx brain tissues have NO counterpart here, and that is
# a property of the panel, not of the gene. `match` is carried on every row so a reader can
# restrict to exact anatomy:
#
#   exact      the same structure (GTEx frontal cortex BA9 ~ BrainSEQ DLPFC).
#   adjacent   a different structure in the same division (putamen / nucleus accumbens are
#              striatal like caudate; BA24 is cortex but not dorsolateral prefrontal). A
#              junction may legitimately behave differently here.
#   secondary  DLPFC-matched loci are also carried in caudate, because BrainSEQ's aging
#              switch modules are caudate-derived and that is where the switch layer is
#              best characterized.
#
# Deliberately absent: cerebellum and cerebellar hemisphere (58 of the 76 concordant events,
# and BrainSEQ sequences no cerebellum), hypothalamus, amygdala, substantia nigra, spinal
# cord. Targets in those tissues are emitted with verdict `no_matched_brainseq_region`
# rather than dropped -- silently skipping them would make an untestable locus look like a
# locus that was never nominated.
_TISSUE_TO_REGION: dict[str, tuple[tuple[str, str], ...]] = {
    "Brain_Hippocampus": (("hippocampus", "exact"),),
    "Brain_Cortex": (("dlpfc", "exact"), ("caudate", "secondary")),
    "Brain_Frontal_Cortex_BA9": (("dlpfc", "exact"), ("caudate", "secondary")),
    "Brain_Anterior_cingulate_cortex_BA24": (("dlpfc", "adjacent"),),
    "Brain_Caudate_basal_ganglia": (("caudate", "exact"),),
    "Brain_Putamen_basal_ganglia": (("caudate", "adjacent"),),
    "Brain_Nucleus_accumbens_basal_ganglia": (("caudate", "adjacent"),),
}

_DEEP_DIVE = stage_out("anchoring", "deep_dive", "deep_dive_events.parquet")

# The competing arm each gene's DISPLAY ITEM actually claims, as a coordinate that must
# appear in the event's competing arm. The anchored junction usually takes part in several
# events, each contrasting it with a different alternative form, and they do not all agree;
# without naming the one the figure asserts, "does SNCA validate?" has no single answer and
# the code would be choosing on the reader's behalf.
#
#   SNCA  89838252 -- the canonical distal first exon (ENST00000336904,
#         chr4:89,838,252-89,838,315) that Fig 4A draws the anchored proximal exon against.
#   CTSH  78937686 -- the competing acceptor the same AD risk allele moves the other way,
#         per deep_dive_events.
#   PRDM2 13788079 -- the early 3' terminus (ENST00000413440's last exon end) that Fig 4A
#         draws the distal terminal exon against. PRDM2 replaced SNCA on the panel
#         2026-09-19; SNCA and CTSH are kept here because the falsification text still
#         refers to their contrasts.
_REFERENCE_ARM: dict[str, int] = {"SNCA": 89838252, "CTSH": 78937686, "PRDM2": 13788079}


def out_dir():
    return ensure_dir(stage_out("mechanism", "junction_coloc_confirm"))


# --------------------------------------------------------------------------- #
# Target construction
# --------------------------------------------------------------------------- #
_TARGET_COLUMNS = ["gene_name", "ens", "trait", "tissue", "region", "match",
                   "anchored_junction", "competing_junctions", "n_competing", "clpp",
                   "risk_allele"]


def load_targets(genes: tuple[str, ...] | None = None) -> pd.DataFrame:
    """One row per (gene, trait, tissue, region) anchored junction.

    ``genes`` restricts to named symbols; ``None`` (the default) takes **every** gene whose
    junction lands in the tissue-matched IsoGraph switch pair. The CLI ran on SNCA and CTSH
    alone until 2026-09-20, which meant the arm reported "No testable targets" as soon as
    those two stopped being the display item, and gave no account of the other 28
    concordant genes.

    A target whose GTEx tissue has no BrainSEQ counterpart is still emitted, with
    ``region`` empty and ``match`` ``"none"``; ``run_region`` never sees it and the report
    counts it as ``no_matched_brainseq_region``. Dropping it here instead would make a
    locus BrainSEQ cannot measure indistinguishable from one that was never nominated --
    and since 58 of the 76 concordant events are cerebellar, that is most of the panel.

    ``competing_junctions`` carries the junctions the same risk allele moves the other way,
    where the colocalization reported any; it may legitimately be empty. SNCA is the case
    in point -- the coloc reports only the anchored junction -- and that does NOT make the
    locus untestable, because a PSI event contrasting the anchored junction with an
    alternative form supplies the competitor by construction (``direct`` mode).
    """
    ev = pd.read_parquet(_DEEP_DIVE)
    ev = ev[ev["junction"].astype(str).str.len().gt(0)]
    if genes is not None:
        ev = ev[ev["gene_name"].isin(genes)]

    rows: list[dict] = []
    for (gene, trait, tissue), grp in ev.groupby(["gene_name", "trait", "tissue"]):
        anchored = grp[grp["junction_in_switch_pair"].fillna(False).astype(bool)]
        competing = grp[~grp["junction_in_switch_pair"].fillna(False).astype(bool)]
        if anchored.empty:
            continue
        mapped = _TISSUE_TO_REGION.get(str(tissue), (("", "none"),))
        for region, match in mapped:
            for _, a in anchored.iterrows():
                comp = [str(c) for c in competing["junction"] if str(c).strip()]
                rows.append({
                    "gene_name": gene, "ens": str(a["ens"]), "trait": trait,
                    "tissue": tissue, "region": region, "match": match,
                    "anchored_junction": str(a["junction"]),
                    "competing_junctions": ";".join(comp),
                    "n_competing": len(comp),
                    "clpp": float(a["clpp"]), "risk_allele": str(a["risk_allele"]),
                })
    # Explicit columns: a gene set absent from the deep dive is a legitimate null, and an
    # empty frame without them would crash the per-region filter downstream.
    return pd.DataFrame(rows, columns=_TARGET_COLUMNS)


def unmatched_targets(targets: pd.DataFrame) -> pd.DataFrame:
    """Targets in a GTEx tissue BrainSEQ does not sequence, as report-ready rows."""
    if targets.empty:
        return targets
    u = targets[targets["match"].eq("none")].copy()
    if u.empty:
        return u
    u["mode"] = "none"
    u["verdict"] = "no_matched_brainseq_region"
    u["verdict_reason"] = (
        "BrainSEQ sequences caudate, DLPFC and hippocampus only; "
        + u["tissue"].astype(str) + " has no counterpart, so the junction cannot be "
        "measured in this panel at all")
    return u


def _match_events(junction: str, gene_meta: pd.DataFrame) -> list[str]:
    """Every PSI event of the gene whose coordinates carry this junction's endpoints.

    Coordinates come from ``event_info`` (parsed upstream by ``_load_splice``), never from
    the table's numeric ``start_*``/``end_*`` columns: those are misaligned in the LIBD PSI
    tables -- CTSH rows on chr15 carry the same (100635746, 100636191) tuple that SNCA rows
    on chr4 do -- so matching on them silently returns nothing for every gene.
    """
    parsed = _parse_junction(junction)
    if parsed is None:
        return []
    chrom, s, e = parsed
    hits = []
    for eid, ev_chrom, ev_info in zip(gene_meta["event_id"], gene_meta["chrom"],
                                      gene_meta["event_info"], strict=True):
        if str(ev_chrom) != str(chrom):
            continue
        cs = _event_coords(ev_info)
        if any(abs(x - s) <= _COORD_TOL for x in cs) and any(
            abs(x - e) <= _COORD_TOL for x in cs
        ):
            hits.append(str(eid))
    return hits


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def _residualize(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    """Least-squares residuals of y on design (design already carries an intercept)."""
    ok = np.isfinite(y)
    out = np.full_like(y, np.nan, dtype=float)
    if ok.sum() < design.shape[1] + 3:
        return out
    beta, *_ = np.linalg.lstsq(design[ok], y[ok], rcond=None)
    out[ok] = y[ok] - design[ok] @ beta
    return out


def _design(sample_table: pd.DataFrame, covariates: list[str]) -> np.ndarray:
    cols = [c for c in covariates if c in sample_table.columns]
    X = pd.get_dummies(sample_table[cols], drop_first=True, dummy_na=False)
    X = X.apply(pd.to_numeric, errors="coerce").astype(float)
    X = X.loc[:, X.std(axis=0, skipna=True) > 0]
    X = X.fillna(X.mean(axis=0))
    return np.column_stack([np.ones(len(X)), X.to_numpy(dtype=float)])


def _rho(a: np.ndarray, b: np.ndarray, min_n: int) -> tuple[float, int]:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < min_n:
        return float("nan"), int(ok.sum())
    x, y = a[ok], b[ok]
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return float("nan"), int(ok.sum())
    return float(stats.spearmanr(x, y).statistic), int(ok.sum())


def _pair_rhos(resid: dict[str, np.ndarray], event_ids: list[str], min_n: int,
               exclude: set[str], rng: np.random.Generator,
               max_pairs: int) -> np.ndarray:
    """Spearman rho for within-gene event pairs, excluding the target events."""
    ids = [e for e in event_ids if e not in exclude]
    pairs = [(i, j) for k, i in enumerate(ids) for j in ids[k + 1:]]
    if len(pairs) > max_pairs:
        idx = rng.choice(len(pairs), size=max_pairs, replace=False)
        pairs = [pairs[i] for i in idx]
    vals = [_rho(resid[i], resid[j], min_n)[0] for i, j in pairs]
    return np.asarray([v for v in vals if np.isfinite(v)], dtype=float)


def _empirical_p(obs: float, null: np.ndarray) -> float:
    """One-sided: how often is the null at least as anti-correlated as observed?

    (k + 1) / (n + 1) so a p-value can never be reported as exactly zero from a finite
    null -- the usual permutation convention.
    """
    if not np.isfinite(obs) or null.size == 0:
        return float("nan")
    return float((np.sum(null <= obs) + 1) / (null.size + 1))


def _bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values, NaNs passed through and excluded from the family."""
    q = np.full(pvals.shape, np.nan, dtype=float)
    ok = np.isfinite(pvals)
    if not ok.any():
        return q
    p = pvals[ok]
    order = np.argsort(p)
    ranked = p[order] * p.size / (np.arange(p.size) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def _bh_adjust(res: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """BH over the paired-mode p-values, as one family across the whole panel.

    The pre-registered rule applied `alpha` to a raw p-value, which was written when the
    arm ran on two genes. Over the full concordant panel that is an uncorrected family, so
    the per-target `verdict` is kept exactly as pre-registered and the corrected call is
    added beside it as `verdict_bh`. `direct` rows are not in this family at all -- their
    verdict is a usage threshold, not a test.
    """
    if res.empty or "p_within_gene" not in res.columns:
        return res
    res = res.copy()
    paired = res["mode"].eq("paired") if "mode" in res.columns else pd.Series(
        False, index=res.index)
    for src, dst in (("p_within_gene", "q_within_gene"),
                     ("p_matched_gene", "q_matched_gene")):
        vals = pd.to_numeric(res[src], errors="coerce").to_numpy(dtype=float)
        vals = np.where(paired.to_numpy(), vals, np.nan)
        res[dst] = _bh(vals)
    rho = pd.to_numeric(res.get("rho"), errors="coerce")
    res["verdict_bh"] = np.where(
        paired & rho.lt(0) & res["q_within_gene"].le(alpha)
        & res["q_matched_gene"].le(alpha), "validated",
        np.where(paired, "not_validated", res["verdict"]))
    return res


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run_region(region: str, targets: pd.DataFrame, *, alpha: float, min_n: int,
               n_background_genes: int, max_pairs: int, seed: int,
               min_usage: float, trait: str = "aging",
               ) -> tuple[pd.DataFrame, pd.DataFrame]:
    from isograph.io.artifacts import load_dataset_bundle

    rng = np.random.default_rng(seed)
    tgt = targets[targets["region"] == region]
    if tgt.empty:
        return pd.DataFrame(), pd.DataFrame()

    bundle = load_dataset_bundle(_bundle_path("brainseq", region, trait))
    sample_table = bundle.sample_table.copy()
    sample_ids = [str(s) for s in sample_table["sample_id"]]
    covariates = _covariates("brainseq", trait)

    meta, wide = _load_splice("brainseq", region, sample_ids, genes=None)
    wide = wide.set_index("sample_id")
    shared = [s for s in sample_ids if s in wide.index]
    wide = wide.loc[shared]
    st = sample_table.set_index("sample_id").loc[shared].reset_index()
    design = _design(st, covariates)
    print(f"[{region}] {len(shared)} samples, {meta.shape[0]} PSI events, "
          f"design {design.shape[1]} cols", flush=True)

    # Background genes matched on the number of quantified events: closure strength scales
    # with how many events a gene has, so an unmatched background would confound
    # "this pair is special" with "this gene has few events".
    counts = meta.groupby("gene")["event_id"].size()
    target_genes = set(tgt["ens"].map(lambda s: str(s).split(".")[0]))

    results, nulls = [], []
    for _, row in tgt.iterrows():
        ens = str(row["ens"]).split(".")[0]
        gm = meta[meta["gene"] == ens]
        if gm.empty:
            results.append({**row.to_dict(), "verdict": "gene_not_quantified",
                            "rho": np.nan, "n_samples": 0})
            continue

        a_ids = _match_events(row["anchored_junction"], gm)
        if not a_ids:
            # The gene IS quantified -- this is a coverage gap in the LIBD event
            # catalogue, not absence of the junction. PRDM2 is the case in point: 19
            # events in DLPFC, none reaching past ~13,787,067, so the distal terminal
            # exon the GTEx sQTL tags (chr1:13,816,570-13,823,159) is simply not one of
            # the catalogued events. Saying how many events the gene does have keeps
            # that distinguishable from "gene not measured".
            results.append({**row.to_dict(), "verdict": "junction_not_measured",
                            "mode": "none", "n_anchored_events": 0,
                            "n_gene_events": int(len(gm)),
                            "verdict_reason": (
                                f"gene is quantified ({len(gm)} PSI events) but none "
                                f"carries the anchored junction; a coverage gap in the "
                                f"LIBD event catalogue, not evidence against it")})
            continue

        c_ids: list[str] = []
        for j in str(row["competing_junctions"]).split(";"):
            if j.strip():
                c_ids.extend(_match_events(j, gm))
        c_ids = [c for c in dict.fromkeys(c_ids) if c not in set(a_ids)]

        resid = {e: _residualize(wide[e].to_numpy(dtype=float), design)
                 for e in gm["event_id"] if e in wide.columns}
        info = dict(zip(gm["event_id"], gm["event_info"], strict=True))
        etype = dict(zip(gm["event_id"], gm["event_type"], strict=True))

        # ---- direct mode: events whose two arms carry the switch ----
        # Preferred whenever the anchored junction is one arm of an event, because PSI is
        # then the switch ratio itself and no closure baseline is needed.
        #
        # EVERY candidate event is emitted, not the "best" one. An earlier version kept
        # only the event with the most samples quantified, which silently picked SNCA's
        # contrast against the 89836743 alternative rather than against the canonical
        # distal first exon (89838252) that Fig 4A actually claims -- a selection the
        # reader could not see. Reporting all of them makes the spread visible and takes
        # the choice out of the code.
        direct = [e for e in a_ids if e in wide.columns]

        def _n_ok(e):
            return int(np.isfinite(wide[e].to_numpy(dtype=float)).sum())

        emitted = False
        for e in sorted(direct, key=lambda x: -_n_ok(x)):
            if _n_ok(e) < min_n:
                continue
            v = wide[e].to_numpy(dtype=float)
            v = v[np.isfinite(v)]
            med = float(np.median(v))
            minor = float(min(med, 1.0 - med))
            q1, q3 = (float(x) for x in np.percentile(v, [25, 75]))
            rec = {
                **row.to_dict(),
                "n_anchored_events": len(a_ids), "n_competing_events": len(c_ids),
                "mode": "direct",
                "event_id": e, "event_type": str(etype.get(e, "")),
                "event_info": str(info.get(e, "")),
                "competing_arm_reported": bool(e in set(c_ids)),
                "is_reference_contrast": _is_reference(row["gene_name"],
                                                       str(info.get(e, ""))),
                "n_samples": int(v.size),
                "median_psi": med, "minor_form_usage": minor,
                "psi_iqr": q3 - q1, "psi_q1": q1, "psi_q3": q3,
                "frac_samples_minor_ge_thresh": float(
                    np.mean(np.minimum(v, 1.0 - v) >= min_usage)),
            }
            # Pre-registered: the switch is confirmed only if BOTH forms are actually
            # used. A form at <5% is what the long-read assay already reported for SNCA
            # (0.29%); reproducing that on far more samples is a failure, not a pass.
            rec["verdict"] = "validated" if minor >= min_usage else "not_validated"
            rec["verdict_reason"] = (
                f"minor-form usage {minor:.4f} "
                f"{'>=' if minor >= min_usage else '<'} threshold {min_usage}")
            results.append(rec)
            emitted = True
        if emitted:
            continue

        # ---- paired mode: anchored and competing measured by different events ----
        # Reached only when no single event contrasts the anchored junction (so the
        # direct branch emitted nothing); `rec` is rebuilt here because the direct branch
        # now builds one per event inside its own loop.
        rec: dict = {**row.to_dict(), "n_anchored_events": len(a_ids),
                     "n_competing_events": len(c_ids)}
        if not c_ids:
            rec.update({"mode": "none",
                        "verdict": "no_competitor_measured",
                        "verdict_reason": "no PSI event contrasts the anchored junction, "
                                          "and no competing junction is measured"})
            results.append(rec)
            continue

        cand = [(a, c, _rho(resid[a], resid[c], min_n)[0])
                for a in a_ids for c in c_ids if a in resid and c in resid]
        cand = [(a, c, r) for a, c, r in cand if np.isfinite(r)]
        if not cand:
            rec.update({"mode": "paired", "verdict": "insufficient_samples"})
            results.append(rec)
            continue
        a_best, c_best, rho = min(cand, key=lambda x: x[2])
        n_used = _rho(resid[a_best], resid[c_best], min_n)[1]

        within = _pair_rhos(resid, list(resid), min_n,
                            exclude=set(a_ids) | set(c_ids), rng=rng,
                            max_pairs=max_pairs)
        n_ev = int(counts.get(ens, 0))
        lo, hi = 0.5 * n_ev, 2.0 * n_ev
        pool = [g for g, c in counts.items() if lo <= c <= hi and g not in target_genes]
        rng.shuffle(pool)
        matched: list[float] = []
        for g in pool[:n_background_genes]:
            ge = [e for e in meta.loc[meta["gene"] == g, "event_id"] if e in wide.columns]
            if len(ge) < 2:
                continue
            gr = {e: _residualize(wide[e].to_numpy(dtype=float), design) for e in ge}
            matched.extend(_pair_rhos(gr, ge, min_n, exclude=set(), rng=rng,
                                      max_pairs=max(4, max_pairs // 8)).tolist())
        matched_arr = np.asarray(matched, dtype=float)

        p_within = _empirical_p(rho, within)
        p_matched = _empirical_p(rho, matched_arr)
        validated = bool(
            np.isfinite(rho) and rho < 0
            and np.isfinite(p_within) and p_within <= alpha
            and np.isfinite(p_matched) and p_matched <= alpha
        )
        rec.update({
            "mode": "paired",
            "anchored_event": a_best, "competing_event": c_best,
            "n_candidate_pairs": len(cand), "rho": rho, "n_samples": n_used,
            "within_gene_null_n": int(within.size),
            "within_gene_null_median": float(np.median(within)) if within.size else np.nan,
            "p_within_gene": p_within,
            "matched_gene_null_n": int(matched_arr.size),
            "matched_gene_null_median": (float(np.median(matched_arr))
                                         if matched_arr.size else np.nan),
            "p_matched_gene": p_matched,
            "verdict": "validated" if validated else "not_validated",
            "verdict_reason": f"rho={rho:.3f}, p_within={p_within}, p_matched={p_matched}",
        })
        results.append(rec)
        nulls.append(pd.DataFrame({
            "region": region, "gene_name": row["gene_name"], "trait": row["trait"],
            "null": np.r_[np.full(within.size, "within_gene"),
                          np.full(matched_arr.size, "matched_gene")],
            "rho": np.r_[within, matched_arr],
        }))

    res = pd.DataFrame(results)
    nul = pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame()
    return res, nul


def _write_report(res: pd.DataFrame, alpha: float, min_usage: float,
                  missing: tuple[str, ...] = ()) -> str:
    lines = ["# Short-read junction confirmation of the anchored switch pairs", ""]
    if res.empty:
        why = (f" No anchored switch junction in `deep_dive_events` for: "
               f"{', '.join(missing)}." if missing else "")
        return "\n".join(lines + ["No testable targets." + why])
    # One row per candidate PSI event, so row counts are NOT gene counts: `direct` mode
    # deliberately emits every event that carries the anchored junction rather than
    # picking a best one. The headline is therefore stated per gene.
    n_gene = int(res["gene_name"].nunique())
    untestable = res["verdict"].eq("no_matched_brainseq_region")
    g_untestable = set(res.loc[untestable, "gene_name"]) - set(
        res.loc[~untestable, "gene_name"])
    tested = res[~untestable]
    g_measured = set(tested.loc[tested["mode"].ne("none"), "gene_name"]) if len(
        tested) else set()
    g_val = set(tested.loc[tested["verdict"].eq("validated"), "gene_name"]) if len(
        tested) else set()
    vc = res["verdict"].value_counts()
    lines += [
        f"**{n_gene}** genes carry a junction in the tissue-matched IsoGraph switch pair. "
        f"**{len(g_untestable)}** sit only in GTEx tissues BrainSEQ does not sequence and "
        f"cannot be tested here at all. Of the **{n_gene - len(g_untestable)}** that reach "
        f"a BrainSEQ region, **{len(g_measured)}** have the junction measured by a PSI "
        f"event and **{len(g_val)}** validate on at least one event.", "",
        f"alpha = {alpha}, minor-form usage threshold = {min_usage}. Rows below are one "
        f"per candidate PSI event ({len(res)} rows), not one per gene.", "",
        "Every gene whose junction lands in the tissue-matched IsoGraph switch pair is",
        "tested, not only the ones a figure panel highlights.", "",
        "| verdict | rows |", "| --- | ---: |",
    ]
    lines += [f"| {k} | {int(v)} |" for k, v in vc.items()]

    if len(tested):
        lines += ["", "## Per gene, in the regions BrainSEQ can measure", "",
                  "`events` counts candidate PSI events, `validated` how many of them pass "
                  "the usage threshold. A gene with a mixed count validates on some "
                  "contrasts and not others, which is why the display-item rule below "
                  "names one.", "",
                  "| gene | trait | region | match | events | validated | call |",
                  "| --- | --- | --- | --- | ---: | ---: | --- |"]
        grp = tested.groupby(["gene_name", "trait", "region", "match"], dropna=False)
        for (gene, trait, region, match), sub in grp:
            n_ev = int(sub["mode"].ne("none").sum())
            n_ok = int((sub["verdict"] == "validated").sum())
            call = ("junction not measured" if n_ev == 0 else
                    "validates" if n_ok == n_ev else
                    "mixed" if n_ok else "does not validate")
            lines.append(f"| {gene} | {trait} | {region} | {match} | {n_ev} | {n_ok} | "
                         f"{call} |")
    if "verdict_bh" in res.columns:
        paired_n = int(res["mode"].eq("paired").sum()) if "mode" in res.columns else 0
        if paired_n:
            n_bh = int((res["verdict_bh"] == "validated").sum())
            lines += ["", f"Across the panel the `paired` family ({paired_n} targets) is "
                          f"BH-corrected: **{n_bh}** validate on `verdict_bh`. The "
                          f"uncorrected `verdict` column is the pre-registered two-gene "
                          f"rule and is kept for the record."]
    lines += [
        "",
        "`direct` rows test a single PSI event whose two arms contrast the anchored",
        "junction against the alternative form, so PSI is the switch ratio itself and there",
        "is no compositional-closure confound; the statistic is minor-form usage,",
        "`min(median PSI, 1 - median PSI)`, which does not depend on PSI orientation.",
        "`paired` rows test PSI anti-correlation against the closure baseline.",
        "`no_matched_brainseq_region` rows are concordant events in a GTEx tissue BrainSEQ",
        "does not sequence (chiefly cerebellum); they are untestable here, not negative.",
        "`match` says how close the BrainSEQ region is to the GTEx tissue: exact, adjacent",
        "(same division, different structure) or secondary (caudate carried alongside DLPFC).",
        "",
        "| gene | trait | tissue | region | match | mode | n | minor-form usage | median PSI | rho | q_within | q_matched | verdict |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]

    def f(v, d=3):
        try:
            v = float(v)
        except (TypeError, ValueError):
            return "--"
        return "--" if not np.isfinite(v) else f"{v:.{d}f}"

    def n_(v):
        # Rows that never reached a test carry NaN here, not 0; int(NaN) raises.
        try:
            v = float(v)
        except (TypeError, ValueError):
            return "--"
        return "--" if not np.isfinite(v) else str(int(v))

    for _, r in res.sort_values(
            ["verdict", "gene_name", "region"], kind="stable").iterrows():
        lines.append(
            f"| {r['gene_name']} | {r['trait']} | "
            f"{str(r.get('tissue', '')).replace('Brain_', '')} | "
            f"{r.get('region') or '--'} | {r.get('match', '--')} | "
            f"{r.get('mode', '--')} | {n_(r.get('n_samples'))} | "
            f"{f(r.get('minor_form_usage'), 4)} | {f(r.get('median_psi'))} | "
            f"{f(r.get('rho'))} | {f(r.get('q_within_gene'))} | "
            f"{f(r.get('q_matched_gene'))} | {r['verdict']} |"
        )

    ev = res[res.get("event_info", pd.Series(dtype=str)).notna()] if "event_info" in res else res.iloc[0:0]
    if len(ev):
        lines += ["", "## Events tested", ""]
        for _, r in ev.iterrows():
            lines.append(f"- **{r['gene_name']}** / {r['region']}: `{r['event_info']}` "
                         f"({r.get('event_type', '?')}), n = {n_(r.get('n_samples'))}")

    # The gene-level call rests ONLY on the event the display item claims, so that a
    # gene cannot be declared confirmed on the strength of some other alternative form.
    # .eq(True) rather than .fillna(False): the column is object dtype once untestable
    # rows are concatenated, and fillna on object dtype is deprecated.
    ref = res[res.get("is_reference_contrast",
                      pd.Series(False, index=res.index)).eq(True)
              ] if "is_reference_contrast" in res else res.iloc[0:0]
    lines += ["", "## Display-item contrast (the row the decision rests on)", ""]
    if ref.empty:
        lines.append("No event contrasts the anchored junction against the arm the figure "
                     "claims; the decision rule cannot be applied.")
    else:
        lines += ["| gene | region | event | n | minor-form usage | verdict |",
                  "| --- | --- | --- | ---: | ---: | --- |"]
        for _, r in ref.drop_duplicates(["gene_name", "region", "event_id"]).iterrows():
            lines.append(f"| {r['gene_name']} | {r['region']} | `{r['event_info']}` | "
                         f"{n_(r.get('n_samples'))} | {f(r.get('minor_form_usage'), 4)} | "
                         f"{r['verdict']} |")

    lines += ["", "## Pre-registered decision rule", "",
              "Applies only to genes with a display-item contrast in `_REFERENCE_ARM`; "
              "every other gene is reported above, not decided.", ""]
    ref_genes = (ref.groupby("gene_name")["verdict"]
                 .apply(lambda s: bool((s == "validated").all())).to_dict()
                 if not ref.empty else {})
    for gene, ok in sorted(ref_genes.items()):
        if ok:
            lines.append(
                f"- **{gene} validates.** Report the short-read result as the orthogonal "
                f"confirmation and cite the long-read failure as an assay limitation; "
                f"{gene} may appear on a main figure.")
        else:
            lines.append(
                f"- **{gene} does not validate.** Per the pre-registered rule it stays off "
                f"any main figure and the set-level result stands alone.")
    if not ref_genes:
        lines.append("Decision rule not applicable -- no display-item contrast measured.")

    if "verdict_reason" in res:
        # Untestable-tissue rows all carry the same sentence; one line covers them.
        tested = res[res["verdict"].ne("no_matched_brainseq_region")]
        lines += ["", "Per-target reasons:", ""]
        for _, r in tested.iterrows():
            if isinstance(r.get("verdict_reason"), str):
                lines.append(f"- {r['gene_name']} / {r['region']}: {r['verdict_reason']}")
        un = res[res["verdict"].eq("no_matched_brainseq_region")]
        if len(un):
            by_t = un.groupby("tissue")["gene_name"].nunique().sort_values(ascending=False)
            lines += ["", f"- {len(un)} targets over {un['gene_name'].nunique()} genes sit "
                          f"in a tissue BrainSEQ does not sequence: "
                      + ", ".join(
                          f"{str(k).replace('Brain_', '')} "
                          f"({int(v)} gene{'' if int(v) == 1 else 's'})"
                          for k, v in by_t.items()) + "."]
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", default=None,
                    help="restrict to these symbols; default is EVERY gene whose junction "
                         "lands in the tissue-matched IsoGraph switch pair")
    ap.add_argument("--regions", nargs="+",
                    default=["hippocampus", "dlpfc", "caudate"])
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-n", type=int, default=30,
                    help="minimum samples with both events quantified")
    ap.add_argument("--n-background-genes", type=int, default=300)
    ap.add_argument("--max-pairs", type=int, default=400,
                    help="cap on within-gene event pairs per null (subsampled with --seed)")
    ap.add_argument("--min-usage", type=float, default=0.05,
                    help="minor-form usage a switch must reach to be called "
                         "usable (the long-read failure was 0.0029)")
    ap.add_argument("--seed", type=int, default=13)
    args = ap.parse_args()

    od = out_dir()
    genes = tuple(args.genes) if args.genes else None
    targets = load_targets(genes)
    targets.to_parquet(od / "targets.parquet", index=False)
    testable = targets[targets["match"].ne("none")]
    print(f"{len(targets)} (gene, trait, region) targets over "
          f"{targets['gene_name'].nunique()} genes; {len(testable)} in a BrainSEQ region, "
          f"{len(targets) - len(testable)} in a tissue BrainSEQ does not sequence",
          flush=True)

    res_all, null_all = [], []
    for region in args.regions:
        res, nul = run_region(region, testable, alpha=args.alpha, min_n=args.min_n,
                              n_background_genes=args.n_background_genes,
                              max_pairs=args.max_pairs, seed=args.seed,
                              min_usage=args.min_usage)
        if not res.empty:
            res_all.append(res)
        if not nul.empty:
            null_all.append(nul)

    unmatched = unmatched_targets(targets)
    if not unmatched.empty:
        res_all.append(unmatched)
    res = pd.concat(res_all, ignore_index=True) if res_all else pd.DataFrame()
    res = _bh_adjust(res)
    res.to_parquet(od / "junction_confirm.parquet", index=False)
    if null_all:
        pd.concat(null_all, ignore_index=True).to_parquet(od / "null_rhos.parquet",
                                                          index=False)
    (od / "JUNCTION_COLOC_CONFIRM.md").write_text(
        _write_report(res, args.alpha, args.min_usage,
                      missing=tuple(g for g in (args.genes or ())
                                    if g not in set(targets["gene_name"]))))
    (od / "params.json").write_text(json.dumps(vars(args), indent=2))
    print(f"wrote {od}", flush=True)
    if not res.empty:
        cols = [c for c in ["gene_name", "trait", "region", "match", "mode", "n_samples",
                                    "minor_form_usage", "median_psi", "rho",
                                    "p_within_gene", "q_within_gene", "p_matched_gene",
                                    "q_matched_gene", "verdict"]
                if c in res.columns]
        print(res[cols].to_string(index=False))
        print()
        print(res["verdict"].value_counts().to_string())


if __name__ == "__main__":
    main()
