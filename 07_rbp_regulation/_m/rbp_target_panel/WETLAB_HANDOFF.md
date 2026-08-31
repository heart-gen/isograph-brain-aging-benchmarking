# RBP perturbation — bench hand-off and priority order

**Status: closed for this submission.** The perturbation experiment is designed, costed and
target-selected, but it is a wet-lab study and is out of scope for the Cell Genomics
manuscript. This file is the hand-off: what to run first, what it costs, and what each
tier buys. The full specification is `WETLAB_PERTURBATION_DESIGN.md` in this directory —
read that before ordering anything. Nothing here supersedes it; this is the priority
ordering it does not state.

Everything below is derived from committed, re-runnable analysis:
`rbp_regulon.py` (SLURM 44484238), `rbp_target_panel.py`, `rbp_pair_assayability.py`
(SLURM 44532443). Regenerate with `bash 07_rbp_regulation/_h/07.build_rbp_target_panel.sh --rbp <RBP>`
and `sbatch 07_rbp_regulation/_h/06.rbp_pair_assayability.sh --rbp <RBP>`.

---

## Why there is a priority order at all

The three arms are not three replicates of one question. They span an evidence axis, and
they answer different things, so a partial execution is still informative **provided the
right subset is run**. The single question the experiment adjudicates is whether the
opportunity-adjusted regulon set (43 module x RBP pairs) or the unadjusted hypergeometric
set (714) is the honest one. Any subset that cannot address that is not worth running.

## Priority 1 — NONO + KHDRBS1, both arms, full design

**This is the minimum executable experiment.** Do not run a single arm.

| | |
|---|---|
| conditions | 2 RBPs x 2 siRNAs + scrambled + mock = 6 |
| replicates | n = 4 independent differentiations |
| assay targets | NONO: ATP6V0B, CAPNS1, DDX24, DSTN, PGK1, VDAC3 (6 genes, 11 responsive pairs); KHDRBS1: KXD1, CDIP1, PTPRN (+PPP2R1A, reported separately) |
| matched controls | 10-12 per arm, drawn per design Section 8 |
| readout | isoform-ratio RT-qPCR across the distinguishing junction; every first-choice pair is `junction` class, so one junction-spanning primer per isoform |

**Why this pair and not NONO alone.** NONO is the strongest adjusted arm; KHDRBS1 is the
promiscuity negative control with the widest *unadjusted* footprint and no module surviving
adjustment. Running NONO alone can only produce a confirmation, which is the weakest
possible outcome — it cannot distinguish a real regulon from a generic consequence of
depleting an abundant RBP. The contrast is the experiment. Critically, KHDRBS1 is not
handicapped: it retains 8 responsive pairs over 4 genes against NONO's 11 over 6, so a flat
KHDRBS1 result is evidence of absence rather than absence of assay.

**Built-in internal control, no extra cost.** `DDX24` and `VDAC3` sit in NONO's
`frontal_cortex_ba9/M005` (adjusted OR 0.85, CI 0.69-1.04 — *not* supported by the
adjustment) and are among NONO's best assayable pairs. If NONO knockdown moves
ATP6V0B/CAPNS1/DSTN/PGK1 but not DDX24/VDAC3, the adjusted OR is tracking something real at
module resolution, measured on pairs of comparable assay quality. Exclude them from the
NONO success criterion.

## Priority 2 — add ELAVL1

Adds breadth (11 modules over 6 regions, median adjusted OR 1.68, none contradicted) and
a second independent positive. Add it if the budget allows a third arm from the start;
adding it later is also fine, since the arms are analysed independently.

**Caveat that costs it rank:** only 2 of its 14 tier-1 genes (PPP2R1A, SYP) carry a
measurable *and* responsive pair, against 6 of 13 for NONO. ELAVL1's panel is the thinnest
of the three despite the strongest nomination breadth, so it is the least efficient arm per
plate. `PPP2R1A` is additionally shared with the KHDRBS1 panel and responsive in both — assay
it, but exclude it from the primary regulon-versus-control comparison in both arms.

## Priority 3 — secondary endpoints and follow-on

- Measurable-but-age-static tier-1 pairs as secondary endpoints. They can still move under
  knockdown; age non-response only means the ratio does not drift over the lifespan.
- CLIP in the same model system, **required** before any claim of *direct* regulation.
  NONO is a paraspeckle protein with roles beyond splicing, so a switch response alone does
  not establish direct binding.
- Region-matched model if a specific nomination fails: modules are region-specific and
  targets are assayed in one cell model.

## Do not start without

1. **Baseline expression confirmed** for every assay target in the chosen line —
   `expression_check_list.tsv` per panel directory. Nine of the prespecified pairs are lost
   at this step in silico already; expect further attrition in a cell line.
2. **Protein-level knockdown confirmation** (western, not only RT-qPCR). A failed knockdown
   invalidates that arm's *null*, which is exactly what the KHDRBS1 arm is for.
3. **Prespecified analysis registered** — design Section 9. With a negative-control arm,
   post-hoc reinterpretation of a weak KHDRBS1 response is otherwise trivial.
4. **GENCODE v47** for any sequencing readout, matching the computational annotation, so
   transcript IDs correspond exactly.

## Model system

SH-SY5Y differentiated with retinoic acid (7 d) is acceptable and faster; NGN2-induced
iPSC cortical neurons (DIV 21+) are preferred, because 7 of the 12 modules carrying tier-1
targets are cortical. A non-neuronal line is not an option: the modules are brain-derived,
and the eCLIP support layer is HepG2/K562, so cell-line binding data already cannot
substitute for a neuronal test.

## What the manuscript claims without this experiment

Predicted candidate regulons, not validated ones. The computational layer establishes
motif-differential enrichment adjusted for motif opportunity, plus eCLIP binding *capacity*
at alternative exons in a non-neural cell line. It does not establish that any nominated RBP
regulates its predicted targets in brain. That limitation is stated explicitly in the
manuscript rather than left implied, and this experiment is the named way to close it.
