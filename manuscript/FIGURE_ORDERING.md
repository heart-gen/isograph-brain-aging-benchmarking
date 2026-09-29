# Manuscript figure & table ordering (IsoGraph, Cell Genomics)

The final figure/table sequence for the IsoGraph manuscript, with the one honest claim each
artifact carries. Synthetic-benchmark figures live under `01_synthetic_benchmark/03_metrics/figures/`, where
they are built; every real-data display item lives in `manuscript/_m/figures/`. North-star: IsoGraph is
a **complementary DTU-without-DGE layer**, not a globally superior method — the ordering moves
from "the method recovers switch modules" (synthetic) → "its modules are trustworthy"
(real-data reproducibility) → "they carry real, genetically-anchored biology invisible to
abundance pipelines" (real-data complementarity) → "and the layer is not an artifact of
cell-type composition" (the confound that decides whether any of it is believed).

> **All numbers below are as regenerated on 2026-08-29** from a full re-run of
> `qtl_anchoring` (17 analyses) + `qtl_anchoring_matched` (26) → `qtl_anchoring_meta` →
> `baseline_comparison`. The previously committed downstream tables disagreed with their
> own committed inputs (the `module_enrichment` tables had been regenerated without
> re-running downstream), and the refresh **changed the Fig 3 headline** — see the note
> under Fig 3. Do not quote pre-2026-08-29 values.

## Proposed re-ordering *(suggestion, 2026-09-20 — subject to change)*

> **This is a proposal, not a decision.** The ordering below comes from the results report
> (`reports/results/results_report.md` §7–§8), which organises the paper by biological
> question rather than by analysis stage. The current ordering in **Main figures** stands
> until this is accepted. Nothing downstream has been renumbered.

The report's sequence is *biological motivation → primary finding → the events are real →
genetics as support → the confound that decides belief*. Against the current set that means
**one promotion and one merge**:

| Proposed | Content | Change from current |
| --- | --- | --- |
| **1** | `figConceptOverview` + `fig1_benchmark_overview` — what a switch is, and that co-switch modules are recoverable where switching carries the structure | unchanged |
| **2** | `figTrustFunnel` — modules are donor-robust and their aging axis transfers across cohorts | unchanged |
| **3** | `figOrthogonalConfirm` (+ `figSwitchConsequence`) — the switches are real molecular events, and they are productive UTR/CDS remodeling rather than decay | **promoted** from S-real-8 / S-real-5 |
| **4** | `figGeneticAnchoring` (+ `figQtlSpecificity` folded in or adjacent) — disease variants resolve to specific isoform pairs; PD heritability is splicing-led | **merged**: the specificity contrast becomes a panel or a neighbour rather than its own main figure |
| **5** | `figCompositionRobustness` — the layer is not simply shifting cell-type proportions, regionally | unchanged |

**Why promote orthogonal confirmation.** It is the answer to the first question a reviewer
asks of any short-read isoform result — *is this a quantification artifact?* — and it is
currently answered in the supplement. It is also one of only two genuinely orthogonal layers
in the paper (the other is the allele-aware junction recount), and the strongest version of
the claim that the modules describe real molecular events.

**Why merge the specificity contrast.** Since the 2026-09-19 re-run it is a **representation**
effect, not a method effect — a matched `wgcna_multiplex` baseline on identical features shows
it more strongly (1.084, p = 0.0044 vs IsoGraph 1.065, p = 0.020) — and two of three
pre-specified sensitivity arms no longer clear 0.05. A demoted claim carrying a main figure of
its own invites the question of why it leads. It is still worth showing; it is no longer worth
a numbered slot.

**Costs of the proposal, stated so the decision is informed.** Promoting S-real-8 puts a
figure with two tempering negatives (the continuous statistic is null at P = 0.154; the
genome-wide margin is 0.649 vs 0.623) into the main set — though the current ordering carries
comparable tempering on Figs 3 and 4 already. Merging Fig 3 into Fig 4 makes a dense figure
denser. An alternative that avoids both is to keep five main figures and simply swap Fig 3 and
S-real-8, leaving the merge for later.

**If accepted**, the following must move with it: figure numbers in `content/03.results.md`
and `content/99.supplement.md`, the S-real-N numbering below, the figure-to-claim table in
`reports/results/results_report.md` §8, and the deck in `reports/results/results_briefing.tex`.

---

## Main figures

| Fig | File | Claim |
| --- | --- | --- |
| **1** | `manuscript/_m/figures/figConceptOverview.{pdf,png}` | **What an isoform switch is, and that IsoGraph recovers co-switch modules where truth is known — including the one benchmark that isolates inference from representation.** (A) one gene, two isoforms, usage shifting continuously with age while total gene abundance stays flat — the event a DGE pipeline cannot see. **Redrawn 2026-09-22:** isoform A goes 0.80 → 0.60 and isoform B 0.20 → 0.40, so the shift is continuous but the two isoforms never trade dominance; the earlier crossover at 0.5 read as if a switch *required* a reversal, and a 0.75 → 0.45 / 0.25 → 0.55 shift would still technically cross. (B) the DGE × DTU quadrant map, naming the DTU-without-DGE cell as the paper's contribution; (C) the method in one row (per-gene switch features → VAE latent → gene–gene graph → co-switch modules); (D) module recovery and (E) switch-gene detection over 7,080 synthetic runs, six core scenarios × six methods; **(F, G) replaced 2026-09-22 — the synthetic ground-truth genetics arm** (`genetic_anchoring` scenario: one planted cis-variant per co-switching module plus an equal number of unlinked null variants, four methods on the *identical* switch+abundance features, 120 datasets × 4 methods). (F) module recovery is the decisive separator: IsoGraph VAE 0.997, multiplex 0.982 and Spearman–Leiden 0.997 against **WGCNA 0.466** (Δ = +0.53, one-sided paired Wilcoxon P = 4.9e-19 on all 120 matched datasets — recomputed by the script from the ledger; the 3.0e-18 / n = 117 in the A1 notes predates the completed ledger). (G) genetic-signal detection is **not** a separator and the panel says so: every method recovers ~all planted variants (recall 0.998–1.000) at a calibrated null-variant rate (FPR 0.042–0.054), because a diluted-but-nonzero module still clears Bonferroni at n = 200. The claim the pair supports is that only the switch-aware networks recover the *module* the genetic signal rides on, never that only IsoGraph finds the genetics; the graded best-R² edge (0.458 vs 0.401) is printed by the script and belongs in the caption, not on the axis. The real-data summary panel that closed this figure from 2026-09-19 to 2026-09-22 now opens **Fig 2** (its long-read rows moved to **S-real-8**). Panels A–C are drawn from explicit coordinates and carry no data. Built by `manuscript/_h/concept_overview_figure.R`; 7.09 × 9.5 in, within the 180 × 247 mm limit. The exhaustive 24-panel grid this replaces is retained as built and is the right content for the supplement alongside S1–S12. |
| **2** | `manuscript/_m/figures/figTrustFunnel.{pdf,png}` | **Re-quoted 2026-09-19; panel A added and the rest re-lettered 2026-09-22.** On real brain data IsoGraph's modules are per-module trustworthy: **93/118 chance-trusted across six regions (79%)**, switch drivers reproduce across split halves (median shared-gene driver-loading ρ 0.72–0.82 over all matched pairs, 94.5–100% positive, which is what panel C plots; 0.70–0.85 and 96–100% over reproducible pairs — **re-checked 2026-09-29** against `module_trust_tables/split_half_pairs.csv` on the Leiden-2.0 fits. The 2026-09-20 "correction" to 0.66–0.88 / 0.72–0.89 was itself wrong: it copied the Leiden-5.0 values from the stale `MODULE_TRUST_SUMMARY.md` (see `04_module_trust/README.md`), and the scrambled-join explanation it gave for 0.72–0.82 is not supported), and frozen module signatures carry the same module-specific aging direction in the other cohort — **33/38 BrainSEQ→GTEx (P = 2.1e-6) and 44/55 GTEx→BrainSEQ (P = 4.4e-6)**, where the matched WGCNA arms are null (28/52, P = 0.34; 3/10, P = 0.95). **(A) is the summary panel moved from Fig 1F:** every claim on a common fraction scale beside the matched WGCNA baseline, both counts printed on every row — modules chance-trusted 93/118 vs 62/72; split-half age sign concordance 22/23 vs 14/14; eigengene transfer 33/38 vs 28/52 and 44/55 vs 3/10. It is deliberately unflattering where the data are: IsoGraph sits *below* WGCNA on the first two rows and above it on the transfer rows. The former long-read rows of that panel are not here; they moved to S-real-8. **(B) is not a superiority claim** — WGCNA's trusted rate is 62/72 (86%), *above* IsoGraph's, and the annotation must show both percentages. (C) driver-loading reproducibility. **(D) changed 2026-09-19:** it carries **within-cohort split-half age concordance** (22 of 23 both-significant pairs sign-concordant; WGCNA 14/14), because the cross-cohort module-pair count is null for both methods (IsoGraph 1/38, perm P = 0.91) and was retired as a claim — BrainSEQ is quantified with Salmon and GTEx with RSEM, so that arm confounds processing with biology (PI, 2026-09-19). The cross-cohort attempt moves to `figCrossCohortReplication` in the supplement. (E) structural class of driver switches. 7.2 × 9.0 in. |
| **3** | `manuscript/_m/figures/figQtlSpecificity.{pdf,png}` | **REFRAMED 2026-09-19 on the switching-filter re-run; two earlier claims are withdrawn.** At set level, phenotype-associated switch modules are spared at splicing QTL relative to eQTL (ratio **1.065**, 95% CI 1.010–1.123, p = 0.020, I² = 0.10, k = 10) — about half the legacy effect. **(i) The effect does NOT localise to GO-invisible modules:** GO-invisible 1.080 (p = 0.093) and GO-visible 1.049 (p = 0.126) are both null, and only the phenotype-associated set clears 0.05. **(ii) The effect is NOT IsoGraph-only:** the matched `wgcna_multiplex` baseline, fed identical switch+abundance features, shows a *stronger* contrast (**1.084**, 1.025–1.146, **p = 0.0044**), while `wgcna_switch_only` is null and below 1 (0.927, p = 0.065). The honest claim the panel now carries is a **representation** effect — both methods consuming switch+abundance features show splicing specificity, the switch-only representation does not — and it cannot be used to attribute the effect to IsoGraph's network inference. The legend must also carry the results that run the other way: two of three pre-specified sensitivity arms no longer clear 0.05 (constraint-adjusted 1.039, p = 0.22; SuSiE dose 1.010, p = 0.56), and per gene GTEx gives 62 splicing-only vs 161 expression-only (P = 2.5e-11) and BrainSEQ 16 switch-only vs 94 abundance-only (P = 1.3e-14). Encoding: fill = significance, opacity = arm (all 17 analyses vs shared tissues); k < 3 cells are not drawn. |
| **4** | `manuscript/_m/figures/figGeneticAnchoring.{pdf,png}` | **Rebuilt 2026-09-29 as manuscript Fig 5, panels in citation order: (a) within-donor allelic test (new; shared with `figAllelicImbalance` through `allelic_panels.R`), (b) the 160 CLPP genes by QTL support (82 sQTL only / 18 both / 60 eQTL only) with the 30 resolved genes filled — replaces the verdict bars, whose 60/70/30 split the text could not be reconciled with, (c) the 30 resolved genes by CLPP, (d) PRDM2, (e) LDSC, whose stars are now the one-sided coefficient P (the statistic the text quotes) instead of enrichment P, read from `tableS27_ldsc_partitioned.csv` because the ledger parquet is zstd. The convergence panel E moved to the supplement (`figColocConvergence`); over the 9 testable cells its P range is 0.20–0.86 (LBD/BrainSEQ has no coloc genes), not the 0.19–1.00 quoted below. PD LDSC is splicing-led only in the joint model: the eQTL annotation alone is also enriched (2.93×, coef P = 0.0061, Bonferroni-surviving). The letters in the text below are the pre-2026-09-29 ones.** **Re-quoted 2026-09-19; panel A replaced with PRDM2 the same day.** Disease variants resolve to isoform switches: **(A) PRDM2 (ALS, cortex)** — the colocalizing event is an **alternative 3′ acceptor at a shared donor**, chr1:13,816,570: the risk allele A at rs2744682 lowers usage of the distal junction chr1:13,816,570–13,823,159, which joins PRDM2's distal terminal exon, and raises an unannotated proximal acceptor at 13,821,649, shifting the gene toward an early-terminating form; the junction is in the tissue-matched IsoGraph switch pair. **The alternative arm is drawn as the early-terminating form with PRDM2-206 (ENST00000413440) as its representative** (resolved 2026-09-20). The legend must state that **PRDM2-203 (ENST00000343137) is the same form**, carrying an additional 56-bp 5′ UTR cassette exon and 422 bp less 3′ UTR: both begin at the internal promoter, both terminate ~28 kb before the event's donor so neither can carry the drawn junction, and their CDS is identical. They are not competing partners, and the panel must not be redrawn as a choice between them — PRDM2-206 is the representative only because it is the form ONT can attribute (mean isoform fraction 0.192 vs 2.8e-09) and the one in the single ONT-confirmed switch-like pair against MANE, while the short-read recount favours PRDM2-203 (median minor-form usage 0.108 vs 0.045) because it contributes 3 transcript-specific introns to PRDM2-206's 1. **Do not describe the anchored event as PRDM2's RIZ1/RIZ2 promoter program**: it is a 3′ terminus choice nested inside that program and independent of it, since the distal junction is carried by both RIZ1-type transcripts and by two RIZ2-type ones. PRDM2 was chosen as the only one of the 30 concordant genes corroborated by every genetic layer at once (coloc.abf PP4_sQTL 0.964 vs eQTL 0.494, splicing coloc in 6 tissues and expression coloc in 0; SMR 7/7 supported with no HEIDI rejection and no eQTL instrument; the only concordant event replicating in an independent BrainSEQ cohort). The panel's exact contrast **is** orthogonally supported: in ONT long-read DLPFC BA9 both transcripts are detected and co-expressed with usage Spearman −0.448 and `switch_like = True`, and in the BrainSEQ allele-aware junction recount the drawn junction is an isoform-specific junction of this pair over 292 donors, with 495 of 498 donors carrying both forms. The legend must still state its costs: CLPP is only 0.056, the module is **GO-visible**, junction/switch polarity r is 0.11, and minor-form usage is 0.045, so the early-3′-end form is real but a minority. The LIBD PSI arm returns `junction_not_measured` because that catalogue stops at ~13,787,067 — a gap in that resource only, not a failed confirmation. **SNCA must not be drawn here as resolved** — its LBD junction is not in the tissue's switch pair, it has no PD nomination, and its LBD SMR splits 12 supported against 11 rejected; per the PI decision of 2026-09-18 it is reported as a **falsification example** in the text. **TMEM175 was rejected as the vignette**: its CLPP of 0.483 is contradicted by coloc.abf (PP4_sQTL ~0 in every trait), by zero sQTL/eQTL coloc calls, and by SMR (6/36 supported, its one eQTL signal HEIDI-rejected); TPP1 (CLPP 0.480) survives corroboration but colocalizes on expression too (PP4_eQTL 0.774), so it cannot carry a splicing-specific claim. (B) the aging switch layer carries partitioned heritability; by the coefficient test the splicing annotation is significant in 2 of 5 traits, with **PD splicing-led** (coef p = 0.0050, surviving Bonferroni and the joint model) and **SCZ expression-led** (sQTL τ p = 0.38 vs eQTL τ p = 0.0059). The panel must be per-trait, not pooled. (C, D) of **160** colocalized genes, 82 are sQTL-only, 60 eQTL-only and 18 both; **30 resolve to a specific switch pair** (17 GO-invisible), and only 6 of the 30 reach CLPP ≥ 0.05. BrainSEQ replication is recorded for 21 events over 7 genes but **only 1 of those is also concordant** (PRDM2), so concordance and independent replication are nearly disjoint and the caption must not merge the two counts. Two further caveats bind any single-locus reading: 58 of the 76 concordant events (76%) are cerebellum or cerebellar hemisphere, tracking GTEx sample size rather than disease-relevant tissue, and `n_switch_pairs` sits at its cap of 15 for the median event, so "concordant" means the junction transcript appears in one of up to 15 reported pairs. C–D are the CLPP layer and must be labelled as such — the signal-level hierarchy tiers **121** nominations separately (representative-intron arm), of which **0** recover a known mechanism. **Signal-level re-quote 2026-09-29, all-introns arm (the one the text uses):** 127 analysis × gene rows at PP4_sQTL ≥ 0.8 = **119 gene × trait cells in 116 genes**, **43 sQTL-preferential** (PP4_eQTL < 0.5); headline estimator coloc.abf 89 / coloc.susie 30; paired contrast 67 splicing-only vs 160 expression-only (McNemar P = 5.9e-10). The "42 cells / 40 genes / 14 preferential" in `MANUSCRIPT_PLAN.md` and `VALIDATION_PLAN.md` is the 2026-09-11 run on the five aging analyses only (1,638 rows), superseded by the resolution-2.0 re-run. (E) colocalizing genes do *not* concentrate in modules in any trait (permutation P = 0.19–1.00, null in all 10 cells); this replaced a retracted SCZ panel whose all-genes denominator was ascertainment. Built by `manuscript/_h/genetic_anchoring_figure.R`. |
| **5** | `manuscript/_m/figures/figCompositionRobustness.{pdf,png}` | **Re-quoted 2026-09-19. The switch layer is not simply shifting cell-type proportions — but only partly.** Composition-adjusted counts survive in the limbic/striatal aging arm (BrainSEQ caudate 45→14, DLPFC 22→28) and in 6 of 8 deconvolved GTEx regions (hippocampus 44→29, caudate 87→50, ACC 928→365), while **the two GTEx cortical regions collapse almost completely** (frontal cortex BA9 797→2, cortex 1,169→2) and the SCZD disease signal is largely composition-confounded (60→11, 18% retained). Aging hippocampus contributes 1 gene and nothing to the argument. Confounder-vs-mediator is unresolvable from these data and the legend must say so. Module member genes are not bags of cell-type markers (0 depleted, 1–2 enriched per region). Built by `manuscript/_h/composition_robustness_figure.R`. |

### Fig 3 — what the switching-filter re-run changed *(2026-09-19)*

Two separate corrections have landed on this figure. The 2026-08-29 refresh removed the
GO-invisible localisation; the 2026-09-19 re-run halved the effect and **reversed the
matched-baseline control**. Current values, IsoGraph, full 17-analysis meta:

| module set | k | ratio | p | legacy (expression filter) |
| --- | --- | --- | --- | --- |
| all_modules | 17 | 1.033 | 0.018 | 1.068, 1.3e-5 |
| pheno_sig_modules | 10 | **1.065** | **0.020** | 1.111, 3.6e-4 |
| go_invisible_modules | 7 | 1.080 | 0.093 (n.s.) | 1.068, 0.077 (n.s.) |
| go_visible_modules | 9 | 1.049 | 0.126 (n.s.) | 1.084, 0.050 |

And the matched baselines, which are the reason the figure's claim changed:

| method | features | pheno_sig ratio | p |
| --- | --- | --- | --- |
| isograph | switch+abundance | 1.065 | 0.020 |
| **wgcna_multiplex** | switch+abundance | **1.084** | **0.0044** |
| wgcna_switch_only | switch-only | 0.927 | 0.065 (n.s.) |

**The control reversed.** On the legacy production no matched baseline cleared 0.05 and the
figure carried "the effect is IsoGraph-only". On the switching filter the multiplex baseline,
fed identical features, exceeds IsoGraph. The surviving claim is a *representation* effect:
splicing specificity tracks switch+abundance features regardless of which method consumes
them, and is absent from the switch-only representation. **No sentence in the manuscript may
attribute this contrast to IsoGraph's network inference.**

On the shared-tissue subset neither multiplex method clears 0.05 (IsoGraph 1.040, p = 0.20,
k = 5; wgcna_multiplex 1.059, p = 0.060, k = 5), so that arm is a consistency check, not the
primary.

## Supplementary figures

Synthetic parameter-resolved and stress panels (existing, `01_synthetic_benchmark/03_metrics/figures/`):

| Fig | File | Claim |
| --- | --- | --- |
| S1 | `figS1_idealized_switching` | Module recovery vs switching fraction × noise SD. |
| S2 | `figS2_noise_stress` | Recovery vs count dispersion × noise SD. |
| S3 | `figS3_feature_interactions` | Recovery vs interaction strength × fraction. |
| S4 | `figS4_nonswitching_background` | Recovery & non-switching specificity vs switching fraction. |
| S5 | `figS5_unequal_abundance` | Recovery & switch detection vs isoform-abundance imbalance. |
| S6 | `figS6_scale_compute` | Runtime & peak RAM vs gene count (scale + scale_realistic). |
| S7 | `figS7_abundance_roles` | Abundance-shift detection & isoform-role composition (abundance_switch_mixed). |
| S8 | `figS8_confound_robustness` | Under composition/batch/depth confounds `isograph_vae_residual` holds recovery (AGENTS.md §3). |
| S9 | `figS9_degradation_fallback` | Degradation-aware reliability falls back to the abundance channel as 3′ bias grows. |
| S10 | `figS10_specificity_null` | Non-switching specificity vs the negative-control null. |
| S11 | `figS11_interpretation_accuracy` | Driver/interpretation accuracy on synthetic ground truth. |
| S12 | `figS12_partition_diagnostics` | Partition diagnostics for the synthetic grid (module-count and size distributions per method). Previously built but unnumbered — an unnumbered figure in the release bundle reads as an oversight, so it is now cited here. |

Real-data supplements:

| Fig | File | Claim |
| --- | --- | --- |
| S-real-1 | `manuscript/_m/figures/figBaselineRates.{pdf,png}` | **Bounds the claim:** per-module phenotype rate is driven by switch features (wgcna_switch_only 0.389 > isograph 0.274 > wgcna_multiplex 0.187 ≈ wgcna_gene 0.180, post-2026-09-09 baseline re-fit, all 17 stores incl. SCZD; the earlier 0.336 / 0.384 are stale) and GO enrichment by abundance — IsoGraph is not globally superior. **The legend must state that the comparison is feature-matched, not pipeline-matched:** IsoGraph residualizes the discovery covariates inside its fit while the WGCNA baselines regress RIN only, so the rate gap cannot be read as inference alone. |
| S-real-2 | `manuscript/_m/figures/figGwasResolution.{pdf,png}` | MAGMA module-GWAS enrichment is size-confounded; at canonical resolution 5.0 the giant-module artifact disappears (0/8 significant modules are giant vs 18/43 at res 2.0 and 79/99 for gene-level WGCNA) yet the schizophrenia signal survives across six regions. Built by `manuscript/_h/gwas_resolution_figure.R`; summary `05_genetic_anchoring/_m/gwas/GWAS_RESOLUTION_SUMMARY.md`. |
| S-real-3 | `manuscript/_m/figures/figGoInvisible.{pdf,png}` | The GO-invisible SCZD switch modules carry functionally-consequential isoform switching comparable to the genome-wide background, and nearly every member carries a real anticorrelated transcript pair — GO-invisibility is GO's gene-level bias, not low module quality. **Note:** this is a claim about module *content*, and it is unaffected by the Fig 3 refresh; it is now the paper's only support for the GO-invisible framing, since the QTL contrast no longer separates GO-invisible from GO-visible. Built by `manuscript/_h/go_invisible_figure.R`. |
| S-real-4 | `manuscript/_m/figures/figSeparation.{pdf,png}` | **Abundance and isoform structure are separable and the separation adds information:** IsoGraph's per-gene abundance and switch axes are largely orthogonal in every one of the 17 analyses, not only the SCZD store (per-analysis median \|r\| 0.11–0.24, 23–46% of genes \|r\|<0.1; BrainSEQ 0.11–0.13 / 40–46%, GTEx 0.15–0.24 / 23–35%; **A**, faceted by analysis from `03_module_characterization/_m/axis_orthogonality_all.parquet`); the de-confounded incremental test finds composition-unique genes in most cohort/regions whose switch channel carries phenotype signal total abundance misses (**C**); e.g. NREP (`ENSG00000134986`) has flat total abundance across diagnosis (p=0.93) but a significant isoform switch (p=7e-5, **B**). Read alongside Fig 5, which reports how many of these survive cell-type adjustment. Built by `manuscript/_h/abundance_structure_figure.R`. |
| S-real-5 | `manuscript/_m/figures/figSwitchConsequence.{pdf,png}` | The coding consequence of the switch axis is **productive UTR/CDS remodeling, not decay**: across 9 structural classes only UTR-remodeled (1.27×, 10/10 regions) and CDS-remodeled (1.04×, 10/10) are enriched under a within-gene permutation null, while NMD routing, biotype switch and coding-status loss are depleted; the signal is indistinguishable between GO-invisible and GO-visible modules. Built by `manuscript/_h/switch_consequence_figure.R`; summary `06_switch_mechanism/_m/SWITCH_CONSEQUENCE_SUMMARY.md`. |
| S-real-6 | `manuscript/_m/figures/figRbpRegulon.{pdf,png}` | **Secondary annotation only (PI, 2026-09-18); stays in the supplement.** Switch modules carry candidate RBP regulons, but the covariate adjustment **reorders** the nominations rather than thinning them and the figure must show both arms, never their intersection: 384 hypergeometric hits and **416** covariate-adjusted hits over 11,040 module × RBP cells, with only **89 in both** and one cell raw-enriched yet adjusted-*depleted* (Spearman ρ = 0.35 between the two rankings among significant cells). The same pattern holds at motif-family resolution (336 / 467 / 99 shared) and in the intronic scope (488 / 188). A subset of nominated factors show measured **binding capacity** at their regulon genes' switched exons over those genes' own constitutive exons (**B**; 31/35 directionally preferential, **21/35** at BH q ≤ 0.05 over unique nominated genes, median switched−constitutive gap **0.020**, MATR3 the strongest at 0.106). **Framing is deliberately narrow:** ENCODE eCLIP is HepG2/K562, so this is capacity at alternative-exon sequence, NOT neuronal occupancy of these regulons; the neuronal-CLIP validation program is an honest null. **Module ids from the legacy production name different gene sets and must not be carried over** — the previously quoted GTEx BA9 "M008" ELAVL4/KHDRBS1 regulons were dropped, not re-quoted. Built by `manuscript/_h/rbp_regulon_figure.R`; summary `07_rbp_regulation/_m/rbp/RBP_REGULON_SUMMARY.md`. |
| S-real-7 | `manuscript/_m/figures/figClinicalConsequence.{pdf,png}` | Clinical consequence of the switch layer: (A) switch genes are more LoF-constrained than genome-wide in every region (median LOEUF 0.72 vs 0.94; Fisher p≈1e-93); (B) switched exons carry lower ClinVar P/LP density than constitutive exons (ratio 0.18/0.21, robust to CDS-only scope) — expected alternative-exon biology; (C) the colocalized splicing-led genes are themselves constrained. Built by `manuscript/_h/clinical_consequence_figure.R`; summary `06_switch_mechanism/_m/CLINICAL_CONSEQUENCE_META.md`. |
| S-real-8 | `manuscript/_m/figures/figOrthogonalConfirm.{pdf,png}` | **Re-quoted 2026-09-19; the result strengthened.** The switches are not short-read quantification artifacts. On ONT long-read DLPFC (Aguzzoli-Heberle 2024 NBT; Bambu; n = 12 — independent platform, lab, cohort and quantifier) the genetically anchored switch pairs are switch-like at **0.675 against an abundance-matched null of 0.243** (empirical P = 5e-4, 206 detected pairs of 583, 30 splicing-led genes). **Note 2026-09-29:** the 583 / 206 are anchored → partner *combinations*, not distinct transcript pairs — a pair whose two members both carry the junction is counted from each side — so they are 377 / 131 distinct pairs; on distinct pairs the rate is 92/131 = 0.702 vs a 0.241 null (usable 69/94 = 0.734 vs 0.278), P = 5e-4 either way. The null is IsoGraph switch pairs from other genes matched on abundance but **not on detection**; against detected-only background pairs it is 0.642 (P = 0.13) and 0.699 (P = 0.47) on distinct pairs, and **26 of the 30 genes confirm individually** (25 at usable abundance). *(Legacy expression-filter values, not to be quoted: 0.453 vs 0.252, 53 pairs, 12 genes.)* The panel must also carry the two tempering results: the *continuous* statistic does not separate for the anchored set (−0.161 vs a −0.117 null, P = 0.154), and globally the switch set is only slightly above its own null (0.649 vs 0.623, +0.026 over 16,016 pairs). The signal is in the anchored subset, not in switch calls at large — a bare negative usage correlation is close to vacuous as switch evidence. **2026-09-22:** the genome-wide row (0.649 vs its 0.623 matched null, from `global_null_summary.json`) is now drawn in panel A as the first of three sets, beside the anchored sets, so the thin global margin and the wide anchored margin sit on one axis; it and the anchored row were the two long-read rows of the former Fig 1F. Panel B's y-range now follows the gene count (the fixed 12-gene limit had clipped the 30-gene switching-filter result). |
| S-real-11 | `manuscript/_m/figures/figCrossCohortReplication.{pdf,png}` | **NEW 2026-09-19; a documented limitation, not a result.** BrainSEQ→GTEx module aging replication, reported because it was attempted and failed, and because the failure is attributable. (A) paired age effects for matched modules of both methods; (B) the observed concordant-module count against its matching permutation null — **null for both methods** (IsoGraph 1/38, P = 0.91; WGCNA 2/52, P = 0.15), as is the pooled Stouffer arm (Z = −3.87, perm p = 0.14). The caption must carry the reason: BrainSEQ transcripts are quantified with Salmon and GTEx with RSEM, and per-gene switch–age effects transfer across **quantifiers** at Pearson 0.007–0.022 (sign concordance ≈ 0.51) while transferring across **brain regions within one quantifier** at 0.198–0.588 (median 0.321). Changing the quantifier destroys more signal than changing the region, so a module pair failing here has failed a test that confounds processing with biology (PI, 2026-09-19). Within-cohort split-half concordance carries the reproducibility claim in Fig 2D, and cross-cohort eigengene projection — which re-estimates the module axis in the target cohort's own features — carries the transfer claim. Built by `manuscript/_h/crosscohort_replication_figure.R`. |
| S-real-9 | `manuscript/_m/figures/figIsaConcordance.{pdf,png}` | **The switch signal is not an artifact of IsoGraph.** Under satuRn (the DTU test used by IsoformSwitchAnalyzeR v2), IsoGraph's switch genes carry elevated DTU evidence in **10 of 12 testable analyses** (rank-biserial 0.06–0.44; GTEx hippocampus P = 2e-99, cortex P = 1e-60), with 2 null (n. accumbens, spinal cord). **The remaining 5 of 17 analyses are vacuous by construction** — IsoGraph called no switch genes there — and are drawn as explicit not-testable rows so the denominator cannot be misread. The continuous rank-based statistic is the primary test: the binary empirical-FDR overlap is a conservative lower bound, because Efron's empirical null assumes a sparse alternative and pervasive brain DTU violates it. Built by `manuscript/_h/isa_concordance_figure.R`. |
| S-real-10 | `manuscript/_m/figures/figColocConvergence.{pdf,png}` | **Colocalizing switch genes do NOT concentrate in particular modules — in any trait.** Extends the SCZ-only Fig 4E content to all five traits (AD, PD, LBD, ALS, SCZ) on both aging partitions, and adds the two comparisons a count panel cannot make. (A) Against a size-matched null that holds module sizes fixed, no (trait, source) cell is more concentrated than chance (permutation P = 0.19–1.00, 10/10 cells — **re-read 2026-09-29 from `global.parquet`: 0.20–0.86 over the 9 testable cells; LBD/BrainSEQ has no coloc genes and no P**). (B) Per-module counts, with modules too small in the tested pool to be testable drawn as open symbols; no module survives FDR (min q = 0.23 over 46 testable rows), and every colocalizing gene sits at its own locus (0/58 rows are multi-gene single-locus). (C) **Why the published SCZ convergence does not survive.** `scz_age_projection` reports coloc genes enriched in MAGMA-anchored modules at 15/31 = 48% vs a 25% background (hypergeometric P = 0.004), using ALL module genes as the denominator. But a gene can only colocalize if it was CLPP-tested, and anchored modules are MAGMA-enriched for the same GWAS that decides which genes enter the test — so the tested pool is already 45% anchored, and against it 48% is null (P = 0.40). The same inflation appears for AD (P = 0.043 -> 0.17). This is ascertainment, not convergence. Built by `manuscript/_h/coloc_convergence_figure.R`. |
| S-real-12 | `manuscript/_m/figures/figAllelicImbalance.{pdf,png}` | **Added 2026-09-29; manuscript Fig S19, detail behind Fig 5a.** The within-donor allelic test (BrainSEQ allele-aware junction recount, beta-binomial GLMM with a random donor intercept) finds 590 / 247 / 98 gate-family pairs at q < 0.05 (185 / 92 / 37 genes) in caudate / hippocampus / DLPFC, and is calibrated: the same model on lead-homozygous donors gives a nominal-P rate of 0.043–0.059 and λ_GC 1.01–1.07, against λ 2.61–2.98 for the heterozygous test (recomputed from `allelic_test.parquet` and asserted equal to the stage summaries). (a) as Fig 5a; (b) QQ, both arms; (c) GLMM β vs score Z (Spearman 0.86–0.87; 160 pairs with undefined score_z not drawn); (d) within- vs between-donor sign agreement, 59–65% over all fitted pairs and 82–87% at q < 0.05; (e) module cis rate vs age rank, caudate ρ −0.69 (perm P = 0.07, 7 modules) and hippocampus +0.33 (P = 0.29, 14), pooled −0.183 (P = 0.47) — descriptive, and the claim is the null: cis control does not concentrate in the more age-associated modules. Tables S22 / S23. Built by `manuscript/_h/allelic_imbalance_figure.R`. |

### Retired / folded display items

- **`figSczConvergence`** — was built and committed but carried no figure number and could
  not be cited. Briefly folded into **Fig 4E**, then **retired outright on 2026-09-03**
  when that panel was replaced: the result it shows was retracted 2026-08-30 (see the
  Fig 4 row). `manuscript/_h/scz_convergence_figure.R` is kept as the reproducible builder
  for the retracted analysis, but **its output must not be cited as evidence**; the
  per-module candidate RBP regulators remain available in Table S18.

### Deliberately NOT built

- **Module-level genetic anchoring.** 14/177 modules anchored at perm p ≤ 0.05 (vs ~9
  expected), and **not** concentrated in GO-invisible modules (5.6% vs 13% GO-visible).
  Splicing anchoring is a pooled-gene property, not a per-module one; plotting it would
  invite exactly the reading the analysis rules out. It is **Table S16** plus one Results
  sentence naming the exemplar anchored programs, and that is the whole of it.

## Supplementary tables

See `manuscript/_m/supp_tables/SUPPLEMENTARY_TABLES.md` (Tables S1–S18) — three-baseline pooled
(S1) and per-region (S2); QTL specificity contrast (S3), matched-baseline contrast (S4) and
raw cis-QTL ORs (S5); GO-invisible disease modules (S6); per-region trust funnel (S7); the
per-gene deep-dive set backing Fig 4 — verdict panel (S8), per-event anchor→switch→consequence
(S9), per-gene RBP regulators (S10), per-gene/exon clinical annotation (S11), and per-gene
known-isoform-biology literature (S12); and the validation layer — cell-type composition
adjustment (S13, backs Fig 5), long-read orthogonal confirmation (S14, backs S-real-8),
satuRn concordance (S15, backs S-real-9), module-level genetic anchoring (S16, replaces a
figure), per-RBP eCLIP binding support (S17, backs S-real-6B) and SCZ convergence with its
candidate regulators (S18, backs S-real-10; it backed Fig 4E before that panel was replaced). The synthetic benchmark summary —
**supplementary**, moved out of the main text
(`01_synthetic_benchmark/03_metrics/_m/tableS_benchmark_summary.csv`; six core accuracy scenarios ×
six main methods) — and the scale-compute table (`tableS_scale_compute_summary.csv`) accompany
Figs 1/S6.

The integrated **genetic-anchoring Results section** is drafted at
`manuscript/GENETIC_ANCHORING_RESULTS.md`. **It predates the 2026-08-29 refresh and its
GO-invisible localisation claims are now wrong** — rewrite before use.

## Status

**Closed in this pass (2026-08-29):**

1. **Stale downstream outputs re-run.** `qtl_anchoring` (17) + `qtl_anchoring_matched` (26)
   → `qtl_anchoring_meta` → `baseline_comparison`. Fig 3, S-real-1, Table 1 and Tables
   S1–S5 all rebuilt on inputs that now agree with their sources.
2. **Four render defects fixed.** Fig 2C annotation was clipped mid-word, losing both
   permutation p-values; Fig 3A's inset legend sat on top of the GO-visible eQTL point;
   Fig 3's bottom legend truncated "WGCNA multiplex" to "WGCNA mult"; Fig 4B's y-axis
   rendered "Partitioned h.. enrichment" because the superscript-two glyph is absent from
   the export font (now drawn with plotmath).
3. **Two orphans resolved.** `figSczConvergence` → Fig 4E (that panel was itself replaced on 2026-09-03; see below); `figS12` numbered.
4. **Fig 1 rebuilt schematic-first** (`figConceptOverview`) — the paper no longer opens
   with a 24-panel boxplot grid, and it now defines its central object before using it.
5. **Three display items built** from analyses that had none: Fig 5
   (`figCompositionRobustness`), S-real-8 (`figOrthogonalConfirm`), S-real-9
   (`figIsaConcordance`), plus the eCLIP panel folded into S-real-6.
6. **Six supplementary tables added** (S13–S18).
7. **Two overclaiming annotations corrected.** Fig 2A now shows both trusted percentages
   (89% vs 88%) rather than implying a difference the counts do not contain; Fig 4A's
   SNCA panel is documented as illustrative against the repository's own long-read result.

**Stale legends found 2026-09-20 — the main-figure rows were re-quoted on the
2026-09-19 re-run, the supplementary rows and this Status section were not.**

Each was checked against the stage report's §6 table, which is re-read from the result
files. None has been rewritten here except the Fig 2 driver-loading value (itself re-checked
and reverted to 0.72–0.82 on 2026-09-29, see the Fig 2 row); the rest are
listed so the defect is recorded rather than carried into prose.

| Where | Says | Current value | Source |
| --- | --- | --- | --- |
| **S-real-1** | phenotype rate `switch_only 0.389 > isograph 0.274 > multiplex 0.187 ≈ gene 0.180` | **`switch_only 0.214 > multiplex 0.189 > gene 0.178 > isograph 0.146`** — IsoGraph moves from second to **last**, which changes what the panel bounds | stage 03 §6 |
| **S-real-2** | res 5.0 removes the giant-module artifact, `0/8` vs `18/43` at res 2.0, `79/99` WGCNA | **Reversed.** res 5.0 is `20/26` (77%) against res 2.0's `24/37` (65%); WGCNA gene `83/110` (75%). Production moved to res 2.0 on this criterion | stage 02 §6, overview §1b |
| **S-real-5** | UTR `1.27×, 10/10 regions`; CDS `1.04×, 10/10` | **UTR `1.247`, 9/9 regions; CDS `1.040`, 8/9, I² = 0** | stage 06 §6 |
| **S-real-9** | elevated DTU evidence in `10 of 12` testable analyses, `5 of 17` vacuous | **11 of 13** testable (rank-biserial 0.04–0.45), **4 of 17** vacuous; BrainSEQ caudate P = 3.2e-78 | `06_switch_mechanism/_m/isa_concordance/*/summary.json` |
| **Status item 7** (below) | Fig 2A trusted percentages `89% vs 88%` | **79% (93/118) vs 86% (62/72)** — and WGCNA is the higher one | stage 04 §6 |
| **Still-open item 2** (below) | Fig 3's claim is `1.111, p = 3.6e-4` | **`1.065`, p = 0.020**, with the matched `wgcna_multiplex` baseline *stronger* at `1.084`, p = 0.0044 | stage 05 §6 |

S-real-2's is the consequential one: it asserts the opposite of the current production
decision, so a reader following it would defend resolution 5.0.

**Still open:**

1. **Graphical abstract.** Fig 1A–C is the natural basis for one, but a Cell Press
   graphical abstract is a separate artwork and has not been drawn. Whether Cell Genomics
   research Articles require one is still unverified (`MANUSCRIPT_PLAN.md` §11).
2. ~~**Prose must be rewritten for the Fig 3 refresh.**~~ **DONE 2026-09-01.** Every
   drafted assertion of the GO-invisible genetic localisation is rewritten; Fig 3's
   claim is now the phenotype-associated contrast (1.111, p=3.6e-4) with the
   matched-baseline null as its control. ~~Fig 4E (`figColocConvergence` replacement)
   remains open~~ — **DONE 2026-09-03**: panel E now carries the size-matched-null test
   (null in 10/10 cells); see the Fig 4 row above.
3. ~~**`CELL_GENOMICS_IMPACT_PROGRESS.md` Item 2 is stale.**~~ **DONE 2026-09-09.** It
   quoted the deprecated flat-background arm (17/39 supported, median gap 0.013) as the
   headline; requoted to the canonical GC-matched arm the committed
   `rbp_binding_support.parquet` holds — 38 RBPs, 31 preferential, 25 supported, median
   gap 0.024 — with an explicit note naming which arm is which. S-real-6 and Table S17
   read the parquet and were never affected.
