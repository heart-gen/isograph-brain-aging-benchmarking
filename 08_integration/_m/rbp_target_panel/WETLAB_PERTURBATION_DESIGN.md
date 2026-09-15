# RBP perturbation experiment — design specification

**Three-arm knockdown test of the IsoGraph co-switch regulon hypothesis.**
Generated 2026-08-26 from `rbp_regulon.parquet` (SLURM 44484238) and
`rbp_target_panel.py` (`--rank-by adjusted`, the default). Every number below is
traceable to `07_rbp_regulation/_m/rbp/` and `07_rbp_regulation/_m/rbp_target_panel/<RBP>/`.

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
714 module x RBP pairs are significant unadjusted, only 43 survive adjustment. **This
experiment is designed to adjudicate that disagreement**, not merely to confirm a regulon.

## 2. Three arms, chosen to span the evidence axis

| arm | RBP | raw hits | regions | adj. CI>1 | contradicted | median adj. OR | adj. q<=0.05 | role |
|---|---|---|---|---|---|---|---|---|
| **A** | **NONO** | 5 | 3 | 4 | 0 | 1.52 | **2** | focal positive |
| **B** | **ELAVL1** | 11 | 6 | 6 | 0 | 1.68 | **1** | broad positive |
| **C** | **KHDRBS1** | 21 | 8 | 5 | **1** | 1.33 | **0** | promiscuity negative control |

The three-way contrast is the point of the design:

- **NONO** — few modules, but every one strong and internally consistent. Its best module
  (`frontal_cortex_ba9/M008`, adjusted OR 2.09, 95% CI 1.60-2.73, adjusted q = 1.9e-4) is
  the single strongest adjusted result in the whole analysis. eCLIP-supported.
  *Prediction: strong, module-restricted switch response.*
- **ELAVL1** (HuR) — broad **and** strong: 11 modules over 6 regions, median adjusted OR
  1.68, none contradicted. eCLIP-supported. *Prediction: response across several modules.*
- **KHDRBS1** (SAM68) — **the negative control.** It has the widest raw footprint of any
  RBP (21 modules, 8/14 regions), which is why it was the original wet-lab candidate. But
  **no KHDRBS1 module survives opportunity adjustment** (best adjusted q = 0.126 pooled,
  0.102 under the permissive within-RBP correction), its median adjusted OR is 1.33 (rank
  32/59 among RBPs with >=5 raw hits), and its module-level directions are inconsistent —
  5 of 21 modules have OR<1, including `anterior_cingulate_cortex_ba24/M011` at OR 0.558
  (0.380-0.819), i.e. significantly **depleted** after adjustment.
  *Prediction: NO coherent module-restricted response.*

KHDRBS1 is a control for a specific artifact: a degenerate, AU-rich, promiscuous binder
whose motif appears in many switch pairs by chance, producing broad unadjusted
significance without any real shared regulation. If KHDRBS1 knockdown moves its predicted
targets as strongly as NONO's, the opportunity adjustment is over-conservative and the
larger unadjusted regulon set should be trusted. If it does not, the adjusted arm is the
correct read and the manuscript claim narrows accordingly. **Both outcomes are publishable
and the design is powered to distinguish them** — this is not a formality.

## 3. Model system

Human neuronal context is required: the modules are brain-derived (GTEx cortical/subcortical
regions and BrainSEQ DLPFC/caudate/hippocampus), and the eCLIP support layer is HepG2/K562,
**not brain** — so cell-line binding data cannot substitute for a neuronal test.

- **Primary:** SH-SY5Y differentiated to a neuronal phenotype (retinoic acid, 7 d), or
  iPSC-derived cortical neurons (NGN2, DIV 21+) if available. NGN2 neurons are preferred
  for the cortical modules; SH-SY5Y is acceptable and faster.
- **Rationale for a cortical model:** 7 of the 12 modules carrying tier-1 targets across the
  three arms are cortical (`frontal_cortex_ba9`, `cortex`, `anterior_cingulate_cortex_ba24`).
- Confirm baseline expression of every assay target in the chosen line **before** committing
  (`expression_check_list.tsv` in each panel directory lists the genes to check).

## 4. Perturbation

siRNA knockdown, 72 h, is the primary modality; each RBP gets **two independent siRNAs** to
control for off-target effects, analysed separately and required to agree in direction.

| arm | target | expected KD |
|---|---|---|
| A | NONO (ENSG00000147140) | >=70% mRNA, confirmed by protein |
| B | ELAVL1 (ENSG00000066044) | >=70% |
| C | KHDRBS1 (ENSG00000121774) | >=70% |

Controls: non-targeting scrambled siRNA (primary comparator), and mock transfection.
Knockdown efficiency by RT-qPCR **and** western blot; a failed knockdown invalidates that
arm's null result, so protein-level confirmation is mandatory before interpreting arm C.

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
non-regulon controls (Section 8). The tier-1 sets below (13 NONO, 14 ELAVL1, 14 KHDRBS1)
are the nomination universe, not the plate list — Section 7 shows that only 9, 5 and 5 of
them respectively carry a pair an assay can actually resolve.

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

#### KHDRBS1

| tier | gene | ENSG | regions | pairs | adj. OR (95% CI) | module | module GO |
|---|---|---|---|---|---|---|---|
| 1 | **SOCS3** | ENSG00000184557 | 5 | 2 | 3.13 (1.38-7.07) | caudate_sczd/M012 | response to stimulus, protein refolding, regul |
| 1 | **IRF1** | ENSG00000125347 | 5 | 11 | 3.13 (1.38-7.07) | caudate_sczd/M012 | response to stimulus, protein refolding, regul |
| 1 | **OSMR-DT** | ENSG00000249740 | 5 | 11 | 3.13 (1.38-7.07) | caudate_sczd/M012 | response to stimulus, protein refolding, regul |
| 1 | **RPL7** | ENSG00000147604 | 4 | 14 | 1.60 (1.09-2.35) | cortex/M004 |  |
| 1 | **SHROOM3** | ENSG00000138771 | 4 | 13 | 1.60 (1.09-2.35) | cortex/M004 |  |
| 1 | **KXD1** | ENSG00000105700 | 4 | 11 | 1.60 (1.09-2.35) | cortex/M004 |  |
| 1 | **ALOX5AP** | ENSG00000132965 | 4 | 3 | 1.57 (1.09-2.27) | amygdala/M006 | immune system process, immune response, regula |
| 1 | **SLC2A5** | ENSG00000142583 | 4 | 8 | 1.57 (1.09-2.27) | amygdala/M006 | immune system process, immune response, regula |
| 1 | **LAIR1** | ENSG00000167613 | 4 | 10 | 1.57 (1.09-2.27) | amygdala/M006 | immune system process, immune response, regula |
| 1 | **ADAM28** | ENSG00000042980 | 3 | 5 | 1.60 (1.03-2.48) | hypothalamus/M006 | immune system process, immune response, positi |
| 1 | **LAT2** | ENSG00000086730 | 3 | 10 | 1.60 (1.03-2.48) | hypothalamus/M006 | immune system process, immune response, positi |
| 1 | **TNFRSF12A** | ENSG00000006327 | 3 | 5 | 1.60 (1.03-2.48) | hypothalamus/M006 | immune system process, immune response, positi |
| 1 | **MAST1** | ENSG00000105613 | 3 | 12 | 1.37 (1.03-1.83) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 1 | **CEP170B** | ENSG00000099814 | 3 | 2 | 1.85 (0.96-3.53) | cortex/M021 |  |
| 2 | **GABBR2** | ENSG00000136928 | 2 | 30 | 1.37 (1.03-1.83) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 2 | **PTPRN** | ENSG00000054356 | 2 | 13 | 1.37 (1.03-1.83) | frontal_cortex_ba9/M008 | synaptic signaling, trans-synaptic signaling,  |
| 2 | **CDIP1** | ENSG00000089486 | 3 | 9 | 1.33 (0.92-1.93) | frontal_cortex_ba9/M016 | protein folding |
| 2 | **TPCN1** | ENSG00000186815 | 1 | 6 | 1.06 (0.84-1.33) | anterior_cingulate_cortex_ba24/M001 | immune system process, immune response, regula |
| 3 | **MLF2** | ENSG00000089693 | 2 | 10 | 1.85 (0.96-3.53) | cortex/M021 |  |
| 3 | **PPP2R1A** | ENSG00000105568 | 2 | 19 | 1.85 (0.96-3.53) | cortex/M021 |  |
| 3 | **UBE2O** | ENSG00000175931 | 2 | 10 | 1.55 (0.92-2.61) | anterior_cingulate_cortex_ba24/M023 | export from cell, modulation of chemical synap |
| 3 | **TECPR1** | ENSG00000205356 | 2 | 13 | 1.55 (0.92-2.61) | anterior_cingulate_cortex_ba24/M023 | export from cell, modulation of chemical synap |
| 3 | **PPP2R5B** | ENSG00000068971 | 2 | 6 | 1.55 (0.92-2.61) | anterior_cingulate_cortex_ba24/M023 | export from cell, modulation of chemical synap |

### A caveat on three NONO tier-1 genes

`DDX24`, `ATP5F1B` and `VDAC3` are reported under `frontal_cortex_ba9/M005`, whose adjusted
OR is **0.85 (0.69-1.04)** — below 1 and not supported by the adjustment. They entered the
panel on hypergeometric eligibility. Treat them as a **within-arm internal control**: if
NONO knockdown moves the M008/M007 targets but not these, that is direct evidence the
adjusted OR is tracking something real at module resolution. Do not count them toward the
NONO success criterion.

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

| | NONO | ELAVL1 | KHDRBS1 |
|---|---|---|---|
| prespecified pairs | 145 | 129 | 236 |
| structurally measurable | **143 (99%)** | **122 (95%)** | **229 (97%)** |
| — of which `junction` (cheapest design) | 116 | 95 | 193 |
| not discriminable by any internal feature | 2 | 7 | 7 |
| no GTEx quantification (BrainSEQ-only region) | 1 | 19 | 15 |
| measurable (structure + expression) | 48 | 30 | 56 |
| **measurable + responsive** | **11** | **6** | **8** |
| genes with ≥1 measurable + responsive pair | **6 / 15** | **2 / 20** | **4 / 24** |

**Structure is not the bottleneck; expression is.** Nearly every prespecified pair is
resolvable in principle — 97% overall. What removes them is that one member is barely
transcribed: of the pairs scored in GTEx but failing, an isoform stays under 1 TPM in *every*
region for 85/94 (NONO), 67/80 (ELAVL1) and 143/164 (KHDRBS1). The motif-differential
partner is frequently a retained-intron or processed-transcript annotation that is real in
the catalogue but near-absent in tissue. This is a property of the switch-pair definition,
not of the RBP nomination, and it applies to all three arms about equally.

### The pairs to assay

One row per gene — the highest-dynamic-range responsive pair. `carrier` names which member
carries the motif; that is the numerator of the ratio. Full per-pair and per-region detail,
including the second- and third-choice pairs, is in `rbp_pair_assayability.tsv` and
`rbp_pair_assayability_by_region.tsv`.

**NONO** (tier-1 throughout)

| tier | gene | transcript 1 / transcript 2 | carrier | assay | best region | median frac | IQR | age β | age q |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **ATP6V0B** | `ENST00000236067.8` / `ENST00000468183.5` | tx2 | junction | hippocampus | 0.18 | 0.18 | +0.033 | 0.0099 |
| 1 | **CAPNS1** | `ENST00000629983.2` / `ENST00000590049.5` | tx2 | junction | frontal_cortex_ba9 | 0.30 | 0.17 | -0.014 | 0.040 |
| 1 | **DDX24** | `ENST00000555054.1` / `ENST00000553400.1` | tx1 | junction | frontal_cortex_ba9 | 0.35 | 0.26 | -0.020 | 0.0035 |
| 1 | **DSTN** | `ENST00000474024.5` / `ENST00000449141.2` | tx1 | junction | frontal_cortex_ba9 | 0.78 | 0.08 | +0.009 | 0.0086 |
| 1 | **PGK1** | `ENST00000476531.1` / `ENST00000491291.1` | tx2 | junction | anterior_cingulate_cortex_ba24 | 0.36 | 0.50 | -0.054 | 0.016 |
| 1 | **VDAC3** | `ENST00000521348.5` / `ENST00000524291.1` | tx1 | junction | frontal_cortex_ba9 | 0.40 | 0.17 | -0.014 | 0.037 |

**ELAVL1**

| tier | gene | transcript 1 / transcript 2 | carrier | assay | best region | median frac | IQR | age β | age q |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **PPP2R1A** | `ENST00000454220.7` / `ENST00000462990.5` | tx2 | junction | frontal_cortex_ba9 | 0.37 | 0.53 | +0.041 | 0.046 |
| 1 | **SYP** | `ENST00000479808.5` / `ENST00000376303.6` | tx1 | junction | anterior_cingulate_cortex_ba24 | 0.66 | 0.39 | +0.049 | 0.00036 |

**KHDRBS1** (negative-control arm)

| tier | gene | transcript 1 / transcript 2 | carrier | assay | best region | median frac | IQR | age β | age q |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **KXD1** | `ENST00000599319.5` / `ENST00000539106.5` | tx2 | junction | frontal_cortex_ba9 | 0.59 | 0.72 | +0.071 | 0.026 |
| 2 | **CDIP1** | `ENST00000562334.5` / `ENST00000399599.7` | tx2 | junction | frontal_cortex_ba9 | 0.62 | 0.56 | +0.072 | 0.0016 |
| 2 | **PTPRN** | `ENST00000295718.7` / `ENST00000443981.5` | tx1 | junction | anterior_cingulate_cortex_ba24 | 0.42 | 0.29 | -0.022 | 0.00064 |
| 3 | **PPP2R1A** | `ENST00000391791.4` / `ENST00000454220.7` | tx1 | junction_and_segment | frontal_cortex_ba9 | 0.64 | 0.53 | +0.047 | 0.016 |

Every first-choice pair but one is `junction` class, so a single junction-spanning primer per
isoform suffices. `PPP2R1A` in the KHDRBS1 arm needs one junction primer and one internal
amplicon.

### Three consequences for the design

**(a) "Assay all tier-1 genes" is not executable as Section 6 states it.** Restricting to
tier 1, the genes with at least one *measurable* pair are 9/13 (NONO), 5/14 (ELAVL1) and
5/14 (KHDRBS1); with a *responsive* pair, 6/13, 2/14 and 1/14. The minimum panel should be
the responsive pairs above, with the remaining measurable-but-static tier-1 pairs assayed as
secondary endpoints — they can still move under knockdown even though age does not move them.

**(b) The negative control is not handicapped.** KHDRBS1 retains 8 responsive pairs over 4
genes, comparable to NONO's 11 over 6. This matters more than it looks: had the KHDRBS1
targets been unmeasurable, a flat KHDRBS1 result would have been uninterpretable — absence of
assay rather than absence of effect. They are measurable and their ratios demonstrably move
with age, so a null under KHDRBS1 knockdown is genuine evidence.

**(c) `PPP2R1A` is shared between the ELAVL1 and KHDRBS1 panels** and is responsive in both.
It cannot contribute to the arm contrast. Assay it, but exclude it from the primary
regulon-versus-control comparison in both arms and report it separately.

One further note on Section 6's caveat: `DDX24` and `VDAC3` — two of the three
adjusted-OR-0.85 internal-control genes — are among NONO's *best* assayable pairs, and
`ATP5F1B` is measurable but static. The internal control is therefore a real test rather
than a technical dead end: if NONO knockdown moves ATP6V0B/CAPNS1/DSTN/PGK1 but not
DDX24/VDAC3, that contrast is measured on pairs of comparable assay quality.

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

Register these before unblinding. The negative-control arm makes prespecification essential
— otherwise a weak KHDRBS1 response is trivially reinterpretable as a positive.

**Per gene:** Δ(isoform ratio) = log2(ratio_knockdown / ratio_control), averaged over the two
siRNAs (required to agree in sign; disagreement = that gene is uninformative, report as such).

**Per arm:** compare mean |Δ| in tier-1 regulon genes versus matched controls by
Mann-Whitney U (genes are the unit, n ~ 13 vs ~ 11). Report effect size and CI, not only p.

**Success criteria:**

| arm | criterion for "regulon confirmed" |
|---|---|
| A — NONO | tier-1 |Δ| significantly > matched controls (p<0.05), **and** the M008 targets (MATK, KIAA0513, GABBR2) shift consistently in sign |
| B — ELAVL1 | tier-1 |Δ| significantly > matched controls (p<0.05) |
| C — KHDRBS1 | **expected to FAIL the above.** Passing it falsifies the opportunity adjustment |

**Interpretation grid:**

| A/B | C | conclusion |
|---|---|---|
| pass | fail | Adjusted arm is correct. Report the 43 adjusted-significant regulons; the 714 unadjusted are opportunity-inflated. **Expected outcome.** |
| pass | pass | Adjustment is over-conservative. Broader unadjusted regulon set is defensible; revisit the covariate model (it may absorb real length/composition-linked biology). |
| fail | fail | No arm supports the trans-regulon hypothesis at this power. Report as a negative result; the module coordination mechanism is not shared-RBP binding, or the model system is wrong. |
| fail | pass | Anomalous — suspect knockdown efficiency, model-system mismatch, or assay artifact before reinterpreting. |

**Power.** With n=4 and ~13 vs ~11 genes per arm, this design detects a difference of
roughly 1 standard deviation in |Δ| between regulon and control sets at ~80% power. It is
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
  bound*; this is precisely why arm C exists as an empirical check rather than an assumption.
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
| target panels | `rbp_target_panel.py --rbp {NONO,ELAVL1,KHDRBS1}`, rank-by adjusted |
| switch pairs / primer targets | `07_rbp_regulation/_m/rbp_target_panel/<RBP>/rbp_target_switch_pairs.parquet` |
| expression pre-check | `07_rbp_regulation/_m/rbp_target_panel/<RBP>/expression_check_list.tsv` |
| pair assayability + age response | `rbp_pair_assayability.py --rbp {NONO,ELAVL1,KHDRBS1}`, SLURM 44532443, 2026-08-26 |
| — isoform ratios | GTEx v11 RSEM transcript TPM, 13 brain regions (`inputs/processed/gtex_v11/`) |
| — age model | logit(carrier fraction) ~ AGE + SEX + SMRIN + SMTSISCH, OLS, BH within RBP |
| annotation | GENCODE v47 — use the same release for RNA-seq quantification |

Regenerate any panel with:

```bash
bash 07_rbp_regulation/_h/07.build_rbp_target_panel.sh --rbp NONO
sbatch 07_rbp_regulation/_h/06.rbp_pair_assayability.sh --rbp NONO
```
