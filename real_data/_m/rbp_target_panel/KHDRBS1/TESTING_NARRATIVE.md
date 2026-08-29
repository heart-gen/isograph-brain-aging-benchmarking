# KHDRBS1 perturbation: testing narrative and target re-ranking

Companion to `RBP_TARGET_PANEL.md`. That file ranks candidates on the **hypergeometric**
regulon q-value. This file re-ranks the same candidates on the **covariate-adjusted**
(length / GC / region-composition / n_transcripts) GLM, states what the experiment can and
cannot conclude, and specifies the controls the design needs to stay interpretable.

Sources: `real_data/_m/rbp/rbp_regulon.parquet`, `rbp_motif_families.parquet`,
`rbp_counts.parquet`, `rbp_binding_support.parquet`, `real_data/_m/deep_dive/deep_dive_rbp.parquet`,
`rbp_target_switch_pairs.tsv`.

---

## 1. What the statistics do and do not license

KHDRBS1's best covariate-adjusted result is **q_glm = 0.0685** (frontal cortex BA9 / M008,
OR **1.37 [1.03, 1.83]**), over 214 module x region tests and a 34k-test BH correction.

For **target prioritization** this is a reasonable basis: the CI excludes 1, the direction is
consistent with the unadjusted enrichment, and holding follow-up selection to a discovery
threshold is the wrong standard — the perturbation *is* the test. Fixing the separated-fit bug
(task #11) moves it to 0.0641, so the number is stable rather than an artifact of the BH
denominator.

For the **manuscript** the wording must stay at "prioritized candidate regulon", not "KHDRBS1
regulon". The claim the data support is: *a set of co-switching modules is enriched for
KHDRBS1 motif-differential switch pairs, and that enrichment is directionally robust to
opportunity adjustment in the cortical synaptic modules.*

---

## 2. Adjusted evidence for the 10 sampled modules

The panel drew 3 genes each from 10 (region, module) regulons. Under adjustment they separate
into three groups:

| region / module | hyper enr | q_hyper | **OR_adj** | 95% CI | verdict |
|---|---|---|---|---|---|
| caudate_sczd / M012 | 1.72 | 8.7e-03 | **3.13** | [1.38, 7.07] | adjusted-supported |
| hypothalamus / M006 | 1.52 | 2.9e-06 | **1.60** | [1.03, 2.48] | adjusted-supported |
| amygdala / M006 | 1.60 | 5.9e-06 | **1.57** | [1.09, 2.27] | adjusted-supported |
| frontal_cortex_ba9 / M008 | 1.82 | 5.0e-20 | **1.37** | [1.03, 1.83] | adjusted-supported |
| anterior_cingulate / M023 | 1.95 | 2.5e-06 | 1.55 | [0.92, 2.61] | consistent, CI spans 1 |
| frontal_cortex_ba9 / M023 | 1.80 | 3.4e-05 | 1.42 | [0.86, 2.35] | consistent, CI spans 1 |
| anterior_cingulate / M013 | 1.50 | 9.2e-04 | 1.12 | [0.77, 1.63] | consistent, CI spans 1 |
| anterior_cingulate / M001 | 1.19 | 4.8e-02 | 1.06 | [0.84, 1.33] | consistent, CI spans 1 |
| cerebellar_hemisphere / M001 | 1.25 | 2.0e-05 | **0.96** | [0.78, 1.18] | **direction reverses** |
| frontal_cortex_ba9 / M005 | 1.45 | 5.3e-11 | **0.92** | [0.74, 1.14] | **direction reverses** |

**Six panel genes should be dropped.** `frontal_cortex_ba9/M005` (GDI1, JPT1, MDH2) and
`cerebellar_hemisphere/M001` (CDIP1, GNB2, PUS1) both have adjusted point estimates *below* 1.
Their strong hypergeometric q-values (5.3e-11, 2.0e-05) come from transcript length and GC
content, not from KHDRBS1 motif preference. A perturbation result in these genes would be
uninterpretable in either direction.

Note this is not a general deflation — the adjustment *raises* caudate_sczd/M012 from 1.72 to
3.13. It reorders the panel rather than flattening it.

---

## 3. Testability by model system

The four adjusted-supported modules split cleanly by biology, and this decides the experiment:

- **Immune / inflammatory** — hypothalamus/M006, amygdala/M006, caudate_sczd/M012.
  Genes: SOCS3, IRF1, OSMR-DT, CSF3, P2RY6, FSTL3, ZFP36. These carry the *strongest* adjusted
  ORs (1.57–3.13) but are **not testable in a neuronal culture**. They require microglial or
  astrocytic models, or a mixed system.
- **Synaptic** — frontal_cortex_ba9/M008 (OR 1.37 [1.03, 1.83]), with
  anterior_cingulate/M023 (OR 1.55) as the concordant second region.
  Genes: MAPRE3, GABBR2, PTPRN, CEP170B.

The gene deep dive converges independently on the same synaptic set: GABBR2, PTPRN and CHRNB2
all appear in both `frontal_cortex_ba9/M008` and `anterior_cingulate_cortex_ba24/M023`. Two
cortical regions nominating the same genes through the same module is the strongest internal
replication available here.

---

## 4. Recommended core panel

If the validation system is neuronal, test the **cortical synaptic regulon**:

| gene | module(s) | adj. OR | motif-differential pairs | distinct motif carriers | ELAVL1-diff? | orthogonal anchor |
|---|---|---|---|---|---|---|
| **GABBR2** | fcx/M008 + ACC/M023 | 1.37 / 1.55 | 30 | 4 | no | SCZ sQTL coloc |
| **PTPRN** | fcx/M008 + ACC/M023 | 1.37 / 1.55 | 13 | 6 | no | ALS sQTL coloc |
| **MAPRE3** | fcx/M008 | 1.37 | 13 | 3 | no | recurs in 3 regions |
| CEP170B | fcx/M023 | 1.42 | — | — | no | GO-invisible module |

All four are KHDRBS1-motif-differential and **not** ELAVL1-motif-differential, which is what
makes the parallel ELAVL1 arm in Section 6 interpretable on these same genes.

**CHRNB2 has 0 motif-differential switch pairs** despite appearing in the deep dive — there is
no assay readout for it. Do not include it.

If a glial or mixed system is available, run the immune arm in parallel — it has the better
statistics (SOCS3/IRF1/OSMR-DT at OR 1.60, ZFP36 at 3.13) and would be the stronger result.

---

## 5. The eCLIP column is cross-tissue and should not drive ranking

`RBP_TARGET_PANEL.md` uses eCLIP support as its second ranking key (1103 of 1187 candidates
"eCLIP-supported"). Three reasons this carries less weight than its position implies:

1. **Not brain.** ENCODE eCLIP is HepG2/K562 (documented at `rbp_binding.py:14`). It measures
   binding *capacity* in a non-neural cell line, not occupancy in cortex. The project's
   neuronal-CLIP suite covers only TARDBP, NOVA1, NOVA2, RBFOX2 and PTBP2 — **there is no
   brain or neuronal CLIP for KHDRBS1 anywhere in this project.**
2. **Near-saturation, so barely discriminating.** KHDRBS1 has 706,146 peaks across 3 files.
   It marks 57.6% of switched genes and 53.9% of constitutive ones — a 3.7-point difference,
   near the bottom of the panel. A filter that admits ~1103/1187 candidates is not a filter.
3. **Length-driven.** Panel-wide, 26 of 38 RBPs show McNemar OR > 1 with a length-normalised
   density ratio CI entirely below 1; KHDRBS1's density ratio is 0.748 [0.708, 0.788]. Since
   this is K562/HepG2 it is weak evidence *against* a brain program — but it is also not
   evidence *for* one, which is how the column currently functions.

**This strengthens the case for the experiment rather than weakening it.** Brain-specific
KHDRBS1 binding and regulatory evidence does not currently exist for these targets; the
perturbation supplies exactly the missing layer. It should be framed that way in the paper —
as generating the orthogonal evidence, not confirming it.

If a binding readout is wanted alongside perturbation, KHDRBS1 CLIP in a neural system would
be the highest-value addition, and would let these targets be scored the way the NOVA2 arm
already is.

---

## 6. Run a parallel ELAVL1 knockdown — it is the design's load-bearing arm

The panel carries 8 KHDRBS1 matrices spread across **4 motif families**: F0114 (88 matrices /
36 RBPs, ELAVL1-labelled), F0054 (38 / 27, PABPC1), F0112 (47 / 19, ELAVL4), F0058 (8 / 4,
PCBP1). Four of the KHDRBS1 matrices sit in F0114 at `max_sim_within = 1.000` — their PWM is
identical to other matrices in a 36-RBP family. A KHDRBS1 knockdown alone therefore cannot
distinguish KHDRBS1 from an ELAVL-family protein acting at the same U/A-rich sites.

### 6a. ELAVL1 is not merely a control — it is the better-supported candidate

In the same modules, on the same adjusted model, ELAVL1 dominates KHDRBS1 everywhere:

| region / module | KHDRBS1 OR_adj [CI] (q_glm) | **ELAVL1 OR_adj [CI] (q_glm)** |
|---|---|---|
| frontal_cortex_ba9 / M008 | 1.37 [1.03, 1.83] (0.374) | **1.82 [1.36, 2.44] (0.016)** |
| anterior_cingulate / M023 | 1.55 [0.92, 2.61] (0.580) | **2.22 [1.31, 3.76] (0.138)** |
| frontal_cortex_ba9 / M023 | 1.42 [0.86, 2.35] (0.690) | **1.69 [1.01, 2.82] (0.443)** |

ELAVL1 is one of only two RBPs in the whole panel to clear covariate-adjusted BH anywhere
(q_glm = 0.016 in fcx/M008). It is also less composition-sensitive (flat/composition ratio
1.47 vs KHDRBS1's 2.28, panel median 1.16) and more discriminating in eCLIP (rate_diff 0.054
vs 0.037). Framing the second arm as a "specificity control" understates it: on the current
evidence ELAVL1 is the stronger hypothesis and KHDRBS1 the weaker one.

### 6b. The chosen targets separate the two proteins

Despite the PWM overlap, the motif-differential calls do separate. Genome-wide, 74.2% of
KHDRBS1-differential gene x region calls (7,745 / 10,437) are **not** ELAVL1-differential. And
all four recommended targets fall in that discriminating set, in both cortical regions:

| gene | fcx/M008 | ACC/M023 | KHDRBS1-differential | ELAVL1-differential |
|---|---|---|---|---|
| GABBR2 | yes | yes | **yes** | no |
| PTPRN | yes | yes | **yes** | no |
| MAPRE3 | yes | yes | **yes** | no |
| CEP170B | fcx/M023 | ACC/M031 | **yes** | no |

So the two-arm knockdown is a genuine discrimination on these targets, not a formality.

### 6c. The 2x2 sits inside a single module

`frontal_cortex_ba9/M008` (273 genes) partitions into all four cells, which means the whole
design runs in one module with co-regulated background held constant:

| cell | n genes | examples | predicted KHDRBS1 KD | predicted ELAVL1 KD |
|---|---|---|---|---|
| KHDRBS1-only | **89** | GABBR2, PTPRN, MAPRE3 | shift | no shift |
| both | **90** | CELSR3, RTN4R, SNCB, ADGRL1, ACTL6B | shift | shift |
| ELAVL1-only | **7** | ENSG00000099849, ENSG00000106327, ENSG00000108309, ENSG00000149654, ENSG00000160963, ENSG00000167700, ENSG00000174871 | no shift | shift |
| neither | **87** | — | no shift | no shift |

The `neither` cell is the most valuable control in the design: 87 genes from the *same
co-switching module*, so it separates "KHDRBS1 regulates these switches" from "knocking down
any abundant RBP perturbs this module". The `ELAVL1-only` cell (n = 7, small — treat as
qualitative) is the positive control confirming the ELAVL1 arm worked at all.

### 6d. What this cannot resolve computationally

KHDRBS1's 8 matrices span 4 families, so a KHDRBS1-only call may be driven by its
PABPC1-family or PCBP1-family matrices rather than by the ELAVL1-identical ones. The current
`rbp_counts.parquet` schema does not retain `matrix_id`, so the calls cannot be decomposed per
matrix without a re-scan. This is precisely the ambiguity the parallel arm resolves
empirically — which is the argument for running it rather than deferring it.

Secondary options, if the parallel arm is not affordable: a **rescue** with KHDRBS1
re-expression, or **motif-mutant reporters** for the specific carrier isoform. Neither
separates KHDRBS1 from ELAVL1 as directly.

Related: the KHDRBS1 motif is strongly composition-sensitive — genome-wide hit counts fall
2.28x under the GC-binned background versus a panel median of 1.16x. The `composition`-scope
counts used by the panel are the right ones; do not fall back to the `_flatbg` tables.

---

## 7. Readout must be isoform-pair-level, not gene-level

No recommended target has a single motif carrier: GABBR2 has 4 distinct carriers over 30
pairs, MAPRE3 3 over 13, PTPRN 6 over 13. Measuring total gene expression will dilute or miss
the effect. The readout is the **specific motif-carrier isoform versus its switch partner**,
per `rbp_target_switch_pairs.tsv` (`motif_carrier` column). The predicted direction is that
knockdown shifts usage away from the motif-carrying isoform.

---

## 8. Pre-registered interpretation

Fix this before unblinding, so the result reads cleanly either way. With both arms run, the
outcome is read off the Section 6c cells rather than from the KHDRBS1 targets alone:

- **KHDRBS1-attributable regulation.** KHDRBS1-only genes shift, ELAVL1-only genes do not shift
  under KHDRBS1 KD, and the `neither` cell stays flat in both arms. This is the strong result
  and licenses "KHDRBS1 regulates these co-switching events".
- **Shared-motif regulation.** Both the KHDRBS1-only and `both` cells shift under *either*
  knockdown. Licenses "a U/A-rich-binding protein regulates these switches" — weaker, but still
  a positive mechanistic result, and it would indicate the motif model cannot separate the two
  proteins at these sites.
- **ELAVL1-driven.** The `both` cell shifts under ELAVL1 KD while the KHDRBS1-only cell does
  not shift under KHDRBS1 KD. Given ELAVL1's stronger adjusted statistics (Section 6a) this is
  a live possibility, and it is a publishable positive finding — the target simply changes.
- **Non-specific perturbation.** The `neither` cell shifts too. The experiment is uninformative
  about motif-driven regulation; report as a technical negative.
- **Uninformative regardless of outcome:** the six dropped genes in Section 2.

A negative result on the synaptic arm does not invalidate the module — it bounds the claim to
"co-switching modules enriched for KHDRBS1 motif turnover" without the regulatory mechanism,
which is still publishable as the resource-layer finding.

**Powering note.** The `ELAVL1-only` cell has n = 7 in fcx/M008, too few to carry a formal
test. If the ELAVL1 arm needs a properly powered positive control, draw it from ELAVL1-only
genes pooled across fcx/M008 + ACC/M023 + fcx/M023 rather than from M008 alone.
