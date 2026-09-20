# RBP perturbation experiment — design specification

**Two-arm knockdown test of the IsoGraph co-switch regulon hypothesis.**
Generated 2026-08-26 from `rbp_regulon.parquet` (SLURM 44484238) and `rbp_target_panel.py`
(`--rank-by adjusted`, the default); **amended 2026-09-20** for the switching-filter re-run —
the KHDRBS1 arm is removed (§2, §2a) and Section 7 is re-quoted from SLURM 46448488. Sections
1, 3–6 and 8–11 are otherwise at their 2026-08-26 numbers and are marked where a re-run number
is known to differ. Everything is traceable to `07_rbp_regulation/_m/rbp/` and
`08_integration/_m/rbp_target_panel/<RBP>/`.

---

## 1. The hypothesis being tested

IsoGraph co-switch modules are groups of genes whose isoform usage shifts together with
age. The trans-regulatory hypothesis is that a module is coordinated because its member
genes share a common RNA-binding protein whose binding site is *gained or lost* between
the two switch-pair isoforms. If that is true, depleting the RBP should shift isoform
ratios **specifically in its predicted regulon genes** and not in matched non-regulon genes.

The computational analysis nominates candidate regulons two ways: an unadjusted
hypergeometric test, and a binomial GLM adjusting for motif *opportunity* (transcript
length, GC, 5'UTR/CDS/3'UTR composition, transcript count). The two disagree sharply —
714 module x RBP pairs were significant unadjusted and only 43 survived adjustment, with 43
shared. **The experiment was designed to adjudicate that disagreement**, not merely to confirm
a regulon.

> **Re-quoted 2026-09-20.** On the switching-filter re-run the two sets are **384 raw and 416
> adjusted, sharing 89** — the adjusted arm is now *larger* than the raw one and is a
> reordering of it rather than a subset, so "714, then 43 after adjustment" must not be reused.
> The disagreement the experiment addresses is therefore no longer nesting but *ordering*: which
> cells each method puts at the top. With the promiscuity arm removed (§2a), this experiment can
> no longer settle that on its own.

## 2. Two arms, and the control that was removed

> **Amended 2026-09-20 (PI): the KHDRBS1 arm is removed from the design.** On the
> switching-filter re-run its panel returns **0 measurable-and-responsive switch pairs**
> (against 47 measurable-but-age-static and 6 not assayable, over 134 pairs in 14 genes).
> Section 7(b) below set the condition for keeping a negative control: it is only
> interpretable if its targets are measurable, because otherwise a flat result is absence of
> assay rather than evidence of absence. That condition no longer holds, so the arm is cut
> rather than run uninterpretably. What this costs the design is stated in §2a — it is not
> nothing, and the matched non-regulon controls of Section 8 now carry the specificity
> argument alone.

| arm | RBP | raw hits | regions | adj. CI>1 | contradicted | median adj. OR | adj. q<=0.05 | role |
|---|---|---|---|---|---|---|---|---|
| **A** | **NONO** | 5 | 3 | 4 | 0 | 1.52 | **2** | focal positive |
| **B** | **ELAVL1** | 11 | 6 | 6 | 0 | 1.68 | **1** | broad positive |
| ~~C~~ | ~~KHDRBS1~~ | ~~21~~ | ~~8~~ | ~~5~~ | ~~1~~ | ~~1.33~~ | ~~0~~ | **removed 2026-09-20 — no assayable responsive pair** |

- **NONO** — few modules, but every one strong and internally consistent. Its best module
  (`frontal_cortex_ba9/M008`, adjusted OR 2.09, 95% CI 1.60-2.73, adjusted q = 1.9e-4) is
  the single strongest adjusted result in the whole analysis. eCLIP-supported.
  *Prediction: strong, module-restricted switch response.*
- **ELAVL1** (HuR) — broad **and** strong: 11 modules over 6 regions, median adjusted OR
  1.68, none contradicted. eCLIP-supported. *Prediction: response across several modules.*

### 2a. What the removal costs, stated plainly

KHDRBS1 was not decoration. It was the arm that made the experiment adjudicate rather than
confirm: a degenerate, AU-rich, promiscuous binder whose motif appears in many switch pairs
by chance, producing broad unadjusted significance without real shared regulation. A flat
KHDRBS1 response would have been evidence that the opportunity adjustment is the correct
read; a strong one would have been evidence that the adjustment is over-conservative and the
larger unadjusted set should be trusted.

Without it, **this experiment can confirm predicted regulons but cannot by itself adjudicate
between the adjusted and unadjusted regulon sets.** Two things partly substitute, and neither
fully:

1. **The matched non-regulon controls (Section 8)** — 10-12 per arm, drawn from the same
   modules and matched on assay class and expression. They test specificity *within* an arm,
   which is the question that matters most for a positive result, but they are not a
   promiscuous-binder control.
2. **The internal control inside the NONO arm** — pairs in `frontal_cortex_ba9/M005`, whose
   adjusted OR (0.85, CI 0.69-1.04) does *not* support the module. If NONO knockdown moves
   its supported modules and not M005, the adjustment is tracking something real at module
   resolution, on pairs of comparable assay quality.

If a promiscuity control is wanted later, it must be re-selected on the current panel against
the stated criterion — a candidate with a wide unadjusted footprint, no module surviving
adjustment, **and at least a handful of measurable, age-responsive pairs**. KHDRBS1 met the
first two and fails the third, which is why it is out rather than kept with a caveat. Its
panel and assayability tables are left in `rbp_target_panel/KHDRBS1/` as the record of that
judgement.

## 3. Model system

Human neuronal context is required: the modules are brain-derived (GTEx cortical/subcortical
regions and BrainSEQ DLPFC/caudate/hippocampus), and the eCLIP support layer is HepG2/K562,
**not brain** — so cell-line binding data cannot substitute for a neuronal test.

- **Primary:** SH-SY5Y differentiated to a neuronal phenotype (retinoic acid, 7 d), or
  iPSC-derived cortical neurons (NGN2, DIV 21+) if available. NGN2 neurons are preferred
  for the cortical modules; SH-SY5Y is acceptable and faster.
- **Rationale for a cortical model:** 7 of the 12 modules carrying tier-1 targets across the
  arms are cortical (`frontal_cortex_ba9`, `cortex`, `anterior_cingulate_cortex_ba24`).
- Confirm baseline expression of every assay target in the chosen line **before** committing
  (`expression_check_list.tsv` in each panel directory lists the genes to check).

## 4. Perturbation

siRNA knockdown, 72 h, is the primary modality; each RBP gets **two independent siRNAs** to
control for off-target effects, analysed separately and required to agree in direction.

| arm | target | expected KD |
|---|---|---|
| A | NONO (ENSG00000147140) | >=70% mRNA, confirmed by protein |
| B | ELAVL1 (ENSG00000066044) | >=70% |

Controls: non-targeting scrambled siRNA (primary comparator), and mock transfection.
Knockdown efficiency by RT-qPCR **and** western blot; a failed knockdown invalidates that
arm's null result, so protein-level confirmation is mandatory before interpreting any null.

**Design:** 3 RBPs x 2 siRNAs + 2 controls = 8 conditions, **n = 4 biological replicates**
(independent differentiations/passages, not technical replicates). Randomise plate position
and process all conditions in one batch per replicate to avoid confounding batch with arm.

## 5. Readout

The measured quantity is an **isoform ratio within a gene**, not gene-level expression.
This is essential: the hypothesis is about switching, and a gene-level change would not
test it. Each target gene has a defined switch pair — two transcripts, one carrying the
RBP motif and one not.

- **Primary readout:** isoform-specific RT-qPCR across the distinguishing junction of each
  switch pair, expressed as the ratio motif-carrier : non-carrier. The `motif_carrier`
  column in `rbp_target_switch_pairs.parquet` names which transcript carries the site.
- **Secondary / confirmatory:** targeted RNA-seq or full RNA-seq on a subset, quantified
  with the same GENCODE v47 annotation used computationally, so transcript IDs match exactly.
- Primer design must target the junction or region that **distinguishes** the two
  transcripts. Do not pick the pair by hand: `rbp_pair_assayability.tsv` states, per pair,
  whether a distinguishing junction exists (`assay_class`) and which member carries the
  motif. Section 7 lists the pairs that pass. Pairs whose members differ only at a
  transcript end are flagged `nested_terminal` and are excluded — a proximal-versus-distal
  amplicon measures extension usage, which is a different quantity.

**Direction of effect.** The motif is gained or lost between isoforms; knockdown removes
the protein, not the site. The prediction is that the ratio shifts *toward* the isoform
whose fate depends on that RBP. Direction is gene-specific and is **not** prespecified per
gene — the prespecified quantity is the *module-level consistency* of the shift (Section 9),
which is what the regulon hypothesis actually predicts.

## 6. Targets

Tier 0 is the RBP itself (knockdown control). Tier 1 = highest-confidence regulon members;
tier 2-3 = secondary. `pairs` = number of annotated switch pairs available for primer design.
Each RBP's full ranked universe is in `rbp_target_candidates.tsv`.

**Assay at minimum: every measurable + responsive pair in Section 7**, plus matched
non-regulon controls (Section 8). The tier-1 sets below are the nomination universe, not the
plate list — Section 7, re-quoted 2026-09-20, shows that only 2 of 4 tier-1 NONO genes and 4
of 7 tier-1 ELAVL1 genes carry a pair an assay can resolve, and one each carries a responsive
one. *(The KHDRBS1 nomination table that stood here is removed with its arm; see §2a.)*

#### NONO

| tier | gene | ENSG | regions | pairs | adj. OR (95% CI) | module | module GO |
|---|---|---|---|---|---|---|---|
| 1 | **MATK** | ENSG00000007264 | 2 | 5 | 2.09 (1.60-2.73) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 1 | **KIAA0513** | ENSG00000135709 | 2 | 6 | 2.09 (1.60-2.73) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 1 | **BRWD1** | ENSG00000185658 | 2 | 14 | 1.52 (1.22-1.89) | frontal_cortex_ba9/M007 | cellular localization, intercellular transport |
| 1 | **PGK1** | ENSG00000102144 | 2 | 8 | 1.52 (1.22-1.89) | frontal_cortex_ba9/M007 | cellular localization, intercellular transport |
| 1 | **DSTN** | ENSG00000125868 | 1 | 2 | 1.52 (1.22-1.89) | frontal_cortex_ba9/M007 | cellular localization, intercellular transport |
| 1 | **CAPNS1** | ENSG00000126247 | 2 | 6 | 1.58 (1.20-2.10) | caudate_basal_ganglia/M006 | phospholipase C-activating G protein-coupled r |
| 1 | **DDX10** | ENSG00000178105 | 2 | 10 | 1.58 (1.20-2.10) | caudate_basal_ganglia/M006 | phospholipase C-activating G protein-coupled r |
| 1 | **ST6GALNAC5** | ENSG00000117069 | 2 | 8 | 1.58 (1.20-2.10) | caudate_basal_ganglia/M006 | phospholipase C-activating G protein-coupled r |
| 1 | **ATP6V0B** | ENSG00000117410 | 2 | 12 | 1.48 (1.04-2.10) | cortex/M004 |  |
| 1 | **NDUFA4** | ENSG00000189043 | 2 | 10 | 1.48 (1.04-2.10) | cortex/M004 |  |
| 1 | **DDX24** | ENSG00000089737 | 1 | 10 | 0.85 (0.69-1.04) | frontal_cortex_ba9/M005 | aerobic respiration, proton motive force-drive |
| 1 | **ATP5F1B** | ENSG00000110955 | 1 | 5 | 0.85 (0.69-1.04) | frontal_cortex_ba9/M005 | aerobic respiration, proton motive force-drive |
| 1 | **VDAC3** | ENSG00000078668 | 1 | 6 | 0.85 (0.69-1.04) | frontal_cortex_ba9/M005 | aerobic respiration, proton motive force-drive |
| 2 | **GABBR2** | ENSG00000136928 | 1 | 33 | 2.09 (1.60-2.73) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 2 | **GGNBP2** | ENSG00000278311 | 1 | 10 | 1.48 (1.04-2.10) | cortex/M004 |  |

#### ELAVL1

| tier | gene | ENSG | regions | pairs | adj. OR (95% CI) | module | module GO |
|---|---|---|---|---|---|---|---|
| 1 | **STX1B** | ENSG00000099365 | 4 | 4 | 1.82 (1.36-2.44) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 1 | **CHRNB2** | ENSG00000160716 | 3 | 5 | 1.82 (1.36-2.44) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 1 | **ABCG4** | ENSG00000172350 | 3 | 8 | 1.82 (1.36-2.44) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 1 | **TMEM130** | ENSG00000166448 | 2 | 5 | 2.22 (1.31-3.76) | anterior_cingulate_cortex_ba24/M023 | export from cell, modulation of chemical synap |
| 1 | **SYP** | ENSG00000102003 | 2 | 8 | 2.22 (1.31-3.76) | anterior_cingulate_cortex_ba24/M023 | export from cell, modulation of chemical synap |
| 1 | **CACNB1** | ENSG00000067191 | 2 | 7 | 2.22 (1.31-3.76) | anterior_cingulate_cortex_ba24/M023 | export from cell, modulation of chemical synap |
| 1 | **SBK1** | ENSG00000188322 | 2 | 1 | 1.86 (1.23-2.81) | hippocampus/M012 |  |
| 1 | **MT3** | ENSG00000087250 | 2 | 10 | 1.86 (1.23-2.81) | hippocampus/M012 |  |
| 1 | **C1orf216** | ENSG00000142686 | 2 | 2 | 1.86 (1.23-2.81) | hippocampus/M012 |  |
| 1 | **PPP2R1A** | ENSG00000105568 | 3 | 10 | 2.33 (1.14-4.78) | cortex/M021 |  |
| 1 | **MLF2** | ENSG00000089693 | 2 | 10 | 2.33 (1.14-4.78) | cortex/M021 |  |
| 1 | **AP1M1** | ENSG00000072958 | 2 | 5 | 2.33 (1.14-4.78) | cortex/M021 |  |
| 1 | **CPLX2** | ENSG00000145920 | 2 | 16 | 1.69 (1.01-2.82) | frontal_cortex_ba9/M023 | modulation of chemical synaptic transmission,  |
| 1 | **NACC1** | ENSG00000160877 | 2 | 4 | 1.69 (1.01-2.82) | frontal_cortex_ba9/M023 | modulation of chemical synaptic transmission,  |
| 2 | **TMED4** | ENSG00000158604 | 1 | 4 | 1.39 (0.90-2.16) | caudate/M003 |  |
| 3 | **RINL** | ENSG00000187994 | 2 | 9 | 1.61 (1.01-2.58) | frontal_cortex_ba9/M019 |  |
| 3 | **B4GALNT4** | ENSG00000182272 | 1 | 4 | 1.61 (1.01-2.58) | frontal_cortex_ba9/M019 |  |
| 3 | **ADCY6** | ENSG00000174233 | 1 | 5 | 1.61 (1.01-2.58) | frontal_cortex_ba9/M019 |  |
| 3 | **TMEM63B** | ENSG00000137216 | 2 | 7 | 1.39 (0.90-2.16) | caudate/M003 |  |
| 3 | **CDC37** | ENSG00000105401 | 1 | 5 | 1.39 (0.90-2.16) | caudate/M003 |  |

## 7. Which switch pairs are actually assayable

Section 6 lists the genes; this section lists the **transcript pairs you can put on a plate**.
The two are not the same. A pair earns a place on the panel by carrying a differential motif,
which says nothing about whether an assay can resolve the two transcripts or whether their
ratio has room to move. Both were tested explicitly
(`rbp_pair_assayability.py`, SLURM 44532443, 2026-08-26):

- **Structurally measurable** — each transcript owns a splice junction the other lacks
  (`junction`), or one owns a junction and the other a unique exonic stretch of ≥80 bp
  (`junction_and_segment`). Either way both isoforms can be amplified specifically.
- **Expression-measurable** — both isoforms reach a median of 1 TPM in GTEx brain and the
  carrier fraction sits between 0.10 and 0.90, so a knockdown has headroom to shift it.
- **Responsive** — logit(carrier fraction) regressed on age with SEX, RIN and ischemic time
  gives an age term surviving BH within the RBP (q ≤ 0.05), **and** the fraction's
  interquartile range is ≥ 0.05. The IQR gate is the bench-relevant half: a significant but
  one-point shift is not plate-readable.

Age is the only perturbation observable in human tissue. It establishes that the ratio *is*
regulatable — not that this RBP is what regulates it. That is what the experiment is for.

### Where the pairs are lost

> **Re-quoted 2026-09-20** from `rbp_pair_assayability.parquet` as regenerated on the
> switching-filter re-run (SLURM 46448488, 2026-09-19). The legacy numbers below this line
> were from SLURM 44532443 and do not describe the current panel: the prespecified pair
> counts fell by roughly a third to a half, and **the responsive targets are different genes**.

| | NONO | ELAVL1 | ~~KHDRBS1~~ |
|---|---|---|---|
| prespecified pairs | 53 | 90 | ~~134~~ |
| structurally measurable | **53 (100%)** | **89 (99%)** | ~~128 (96%)~~ |
| — of which `junction` (cheapest design) | 38 | 66 | ~~94~~ |
| not discriminable by any internal feature | 0 | 1 | ~~6~~ |
| no GTEx quantification (BrainSEQ-only region) | 4 | 15 | ~~—~~ |
| measurable (structure + expression) | 20 | 29 | ~~47~~ |
| **measurable + responsive** | **6** | **5** | **~~0~~** |
| genes with >=1 measurable + responsive pair | **3 / 7** | **1 / 9** | **~~0 / 14~~** |

**Structure is not the bottleneck; expression is.** Essentially every prespecified pair is
resolvable in principle. What removes them is that one member is barely transcribed — the
motif-differential partner is frequently a retained-intron or processed-transcript annotation
that is real in the catalogue but near-absent in tissue. This is a property of the switch-pair
definition rather than of the RBP nomination, and it applied to all three arms about equally;
it is why the removed KHDRBS1 arm ended with no responsive pair at all.

### The pairs to assay

One row per gene — the highest-dynamic-range responsive pair. `carrier` names which member
carries the motif; that is the numerator of the ratio. Full per-pair and per-region detail,
including second- and third-choice pairs, is in `rbp_pair_assayability.tsv` and
`rbp_pair_assayability_by_region.tsv`.

**NONO** — 6 responsive pairs over 3 genes

| gene | transcript 1 / transcript 2 | carrier | assay | best region | median frac | IQR | age β | age q |
|---|---|---|---|---|---|---|---|---|
| **CAMK1** | `ENST00000482803.1` / `ENST00000411972.1` | tx1 | junction | anterior_cingulate_cortex_ba24 | 0.62 | 0.34 | +0.038 | 3.6e-4 |
| **VIRMA** | `ENST00000297591.10` / `ENST00000522196.1` | tx1 | junction_and_segment | cerebellum | 0.42 | 0.15 | -0.008 | 0.050 |
| **RNH1** | `ENST00000397604.7` / `ENST00000525522.5` | tx2 | junction | cortex | 0.15 | 0.068 | -0.007 | 0.026 |

**ELAVL1** — 5 responsive pairs, all in one gene

| gene | transcript 1 / transcript 2 | carrier | assay | best region | median frac | IQR | age β | age q |
|---|---|---|---|---|---|---|---|---|
| **NAA10** | `ENST00000464845.6` / `ENST00000393710.7` | tx1 | junction | cerebellum | 0.71 | 0.53 | +0.068 | 0.013 |

### Three consequences for the design

**(a) "Assay all tier-1 genes" is not executable as Section 6 states it, and the gap widened.**
Restricting to tier 1, the genes with at least one *measurable* pair are 2/4 (NONO) and 4/7
(ELAVL1); with a *responsive* pair, 1 and 1. The minimum panel is the responsive pairs above,
with the measurable-but-age-static pairs (20 NONO, 29 ELAVL1 measurable in total) as secondary
endpoints — they can still move under knockdown even though age does not move them.

**(b) The experiment is now narrower than a three-arm adjudication.** ELAVL1's responsive
panel collapses to a single gene, so its arm tests one ratio with several pairs rather than
breadth across modules; NONO carries the design. Read §2a before treating a positive result as
evidence about the regulon set as a whole.

**(c) `PPP2R1A` is no longer a shared-target problem.** It was responsive in both the ELAVL1
and the removed KHDRBS1 panels and had to be excluded from the arm contrast; with KHDRBS1 gone
and PPP2R1A no longer among ELAVL1's responsive pairs, that carve-out lapses.

One further note on Section 6's caveat: the adjusted-OR-0.85 internal-control genes in
`frontal_cortex_ba9/M005` remain the within-arm control for NONO. Verify against the current
`rbp_pair_assayability.tsv` which of them still carries a measurable pair before writing them
into the success criterion — the legacy text named `DDX24` and `VDAC3`, neither of which is in
the current responsive set.

## 8. Matched non-regulon controls

The comparison that carries the result is regulon versus **matched non-regulon**, not
regulon versus zero. For each arm select 10-12 control genes that:

1. are expressed in the model system at comparable level,
2. have an annotated switch pair in the same regions,
3. carry the RBP motif in **neither** isoform of the pair (no site switch), and
4. are matched on transcript length and GC to the tier-1 set (these are the covariates the
   GLM adjusts for — matching on them is what makes the wet-lab test a fair analogue of the
   adjusted analysis).

Draw these from `rbp_target_candidates.tsv` rows with `eligible == False` and no motif
switch. Without this matched set the experiment cannot distinguish a regulon effect from a
generic effect of perturbing any abundant RBP.

## 9. Prespecified analysis and success criteria

Register these before unblinding. With the negative-control arm removed (§2a),
prespecification matters *more* rather than less: the surviving arms can only produce
confirmations, so the comparison against matched non-regulon controls is the only thing
standing between a positive result and a generic consequence of depleting an abundant RBP.

**Per gene:** Δ(isoform ratio) = log2(ratio_knockdown / ratio_control), averaged over the two
siRNAs (required to agree in sign; disagreement = that gene is uninformative, report as such).

**Per arm:** compare mean |Δ| in tier-1 regulon genes versus matched controls by
Mann-Whitney U (genes are the unit; on the current panel n is 4 vs ~11 for NONO and 7 vs ~11
for ELAVL1, which is thin — power should be recomputed before registering). Report effect size and CI, not only p.

**Success criteria:**

| arm | criterion for "regulon confirmed" |
|---|---|
| A — NONO | tier-1 |Δ| significantly > matched controls (p<0.05), **and** the M008 targets (MATK, KIAA0513, GABBR2) shift consistently in sign |
| B — ELAVL1 | tier-1 |Δ| significantly > matched controls (p<0.05) |

**Interpretation grid.** With arm C removed, the grid collapses to one axis and the design
no longer adjudicates the adjusted-versus-unadjusted question by itself (§2a).

| A/B vs matched controls | conclusion |
|---|---|
| pass | The predicted regulon targets respond where matched non-regulon pairs from the same modules do not. Supports the nominated regulon; **does not** by itself decide whether the adjusted or the unadjusted regulon set is the right one, because the promiscuity control is gone. |
| fail | No support for the trans-regulon hypothesis at this power. Report as a negative; either module coordination is not shared-RBP binding, or the model system is wrong. Check knockdown at protein level before concluding. |

**Power.** The legacy calculation assumed ~13 vs ~11 genes per arm at n=4, which detects a
difference of roughly 1 standard deviation in |Δ| at ~80% power. **The current panel is
smaller** — 4 tier-1 NONO genes and 7 ELAVL1, one responsive gene each — so this must be
recomputed before registering; on these counts the arm-level Mann-Whitney is underpowered and
the informative comparison may have to be per pair within gene. It is
**not** powered to detect per-gene effects — do not interpret single genes. If per-gene
resolution is needed, increase to n=6 and pre-register the specific genes.

## 10. Known limitations to carry into interpretation

- **Motif presence is predicted, not measured.** Every target rests on a PWM scan (ATtRACT,
  p<1e-4, GC-binned composition background), not on measured binding in neurons.
- **eCLIP support is HepG2/K562, not brain.** The 25/38 binding-supported figure is binding
  *capacity* at alternative exons (median switched-constitutive gap 0.024, small and
  near-universal), not neuronal occupancy. It cannot corroborate a brain regulon.
- **The covariates correlate with the biology.** An RBP genuinely acting on long, AU-rich
  3'UTRs is partly adjusted away by construction. The adjusted set is therefore a *lower
  bound*. Arm C was the empirical check on this and is removed (§2a), so the point now stands as an assumption the experiment cannot test.
- **Modules are region-specific.** Targets are assayed in one cell model but nominated from
  brain regions; a null in a specific gene may reflect region mismatch rather than a false
  nomination.
- **NONO is a paraspeckle protein** with roles beyond splicing; a switch response does not
  by itself establish a direct binding mechanism. Pair a positive result with CLIP in the
  same model system before claiming direct regulation.

## 11. Provenance

| item | source |
|---|---|
| module x RBP regulons (both arms) | `07_rbp_regulation/_m/rbp/rbp_regulon.parquet`, SLURM 44484238, 2026-08-26 |
| motif scan | `rbp_scan.py`, ATtRACT PWMs, p<1e-4, 20 GC x 2 purine composition bins |
| adjusted GLM | `rbp_regulon.py::_regulon_glm`, covariates log_length, gc, frac_5utr, frac_cds, frac_3utr, log_n_transcripts; BH over 31,452 estimable cells |
| target panels | `rbp_target_panel.py --rbp {NONO,ELAVL1}`, rank-by adjusted |
| switch pairs / primer targets | `07_rbp_regulation/_m/rbp_target_panel/<RBP>/rbp_target_switch_pairs.parquet` |
| expression pre-check | `07_rbp_regulation/_m/rbp_target_panel/<RBP>/expression_check_list.tsv` |
| pair assayability + age response | `rbp_pair_assayability.py --rbp {NONO,ELAVL1}`, SLURM 46448488, 2026-09-19 (legacy: 44532443, 2026-08-26) |
| — isoform ratios | GTEx v11 RSEM transcript TPM, 13 brain regions (`inputs/processed/gtex_v11/`) |
| — age model | logit(carrier fraction) ~ AGE + SEX + SMRIN + SMTSISCH, OLS, BH within RBP |
| annotation | GENCODE v47 — use the same release for RNA-seq quantification |

Regenerate any panel with:

```bash
bash 07_rbp_regulation/_h/07.build_rbp_target_panel.sh --rbp NONO
sbatch 07_rbp_regulation/_h/06.rbp_pair_assayability.sh --rbp NONO
```
