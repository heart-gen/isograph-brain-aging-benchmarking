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

## Main figures

| Fig | File | Claim |
| --- | --- | --- |
| **1** | `manuscript/_m/figures/figConceptOverview.{pdf,png}` | **What an isoform switch is, and that IsoGraph recovers them where truth is known.** (A) one gene, two isoforms, usage crossing over with age while total gene abundance stays flat — the event a DGE pipeline cannot see; (B) the DGE x DTU quadrant map, naming the DTU-without-DGE cell as the paper's contribution; (C) the method in one row (per-gene switch features → VAE latent → gene-gene graph → co-switch modules); (D) module recovery and (E) switch-gene detection over 7,080 synthetic runs, six core scenarios × six methods. Panels A–C are drawn from explicit coordinates and carry no data — they illustrate a definition. Built by `manuscript/_h/concept_overview_figure.R`. The exhaustive 24-panel grid this replaces (`01_synthetic_benchmark/03_metrics/figures/fig1_benchmark_overview`) is retained as built and is the right content for the supplement alongside S1–S12. |
| **2** | `manuscript/_m/figures/figTrustFunnel.{pdf,png}` | On real brain data IsoGraph's modules are per-module trustworthy: 236/266 chance-trusted across six regions (89%), switch drivers reproduce (ρ≈0.77–0.82), and 25/130 modules are concordant for aging cross-cohort (perm P = 0.014). **Panel A is a granularity claim, not a superiority claim** — WGCNA's trusted rate is 64/73 (88%), the same rate at ~3.6× coarser granularity, and the annotation now shows both percentages so the panel cannot be misread. The WGCNA concordance arm is also above its own null (3/53, perm P = 0.040). |
| **3** | `manuscript/_m/figures/figQtlSpecificity.{pdf,png}` | IsoGraph's phenotype-associated switch modules are genetically anchored — splicing QTL are spared relative to eQTL (ratio 1.111, 95% CI 1.048–1.177, p = 3.6e-4, I²=0.23; 1.108 on the 8-tissue matched-baseline set, p = 0.001, I²=0.00) — and **the effect is IsoGraph-only**: the matched WGCNA baselines, which consume identical switch features, are null in every module set (p ≥ 0.41). |
| **4** | `manuscript/_m/figures/figGeneticAnchoring.{pdf,png}` | Disease variants resolve to isoform switches: (A) SNCA risk alleles for LBD and PD both raise usage of the same alternative-first-exon junction, mapping onto one GO-invisible IsoGraph switch pair; (B) the aging switch layer carries partitioned heritability across five traits (splicing- vs expression-lean by trait); (C) 12 splicing-led colocalized genes; (D) of 68 colocalized genes, 12 are splicing-led, 23 splicing-unresolved, 33 expression-led; (E) **RETRACTED 2026-08-30 — this panel must be removed or replaced.** It showed 15/32 SCZ colocalized loci in age-sensitive switch modules against a 25% background (P = 0.0058), but that background is ALL module genes; a gene can only colocalize if it was coloc-tested, and anchored modules are defined by MAGMA SCZ enrichment, so they enter the tested pool preferentially. Against the tested pool the same counts give P = 0.40–1.0. The honest all-trait version, with a size-matched null and both denominators drawn, is **S-real-10** (`figColocConvergence`). Built by `manuscript/_h/genetic_anchoring_figure.R`; summary `05_genetic_anchoring/_m/deep_dive/DEEP_DIVE_SUMMARY.md`. |
| **5** | `manuscript/_m/figures/figCompositionRobustness.{pdf,png}` | **The switch layer is not simply shifting cell-type proportions — but only partly.** Composition-adjusted counts survive in the limbic/striatal aging arm (BrainSEQ caudate 43→17, DLPFC 8→15; GTEx 6/8 deconvolved regions retain) and **collapse entirely in the two GTEx cortical regions** (531→0, 438→0), while the SCZD disease signal is largely composition-confounded (34→2). Confounder-vs-mediator is unresolvable from these data and the legend says so. Module member genes are not bags of cell-type markers (0/3,716 depleted; 33/3,716 enriched). Built by `manuscript/_h/composition_robustness_figure.R`. |

### Fig 3 — what the 2026-08-29 refresh changed

The refreshed contrast **no longer localises the splicing-specificity effect to the
GO-invisible module set.** Committed-vs-refreshed, IsoGraph:

| module set | k (old → new) | ratio (old → new) | p (old → new) |
| --- | --- | --- | --- |
| all_modules | 17 → 17 | 1.068 → 1.068 | 1.4e-5 → 1.3e-5 |
| pheno_sig_modules | 12 → 12 | 1.163 → **1.111** | 3.6e-7 → **3.6e-4** |
| go_invisible_modules | 10 → **11** | 1.172 → **1.068** | 2.3e-5 → **0.077 (n.s.)** |
| go_visible_modules | 11 → 11 | 1.104 → 1.084 | 0.022 → 0.050 |

`all_modules` is bit-identical, which is the control showing the pipeline itself is
unchanged; the GO-invisible arm gained a tissue (k 10 → 11) because the refreshed
`module_enrichment` tables changed GO-visible/GO-invisible set membership. The matched
WGCNA baselines are unchanged and still null.

**What survives:** the genetic-anchoring result and its cleanest internal control — the
effect is real in the phenotype-associated modules and is IsoGraph-only against matched
baselines that consume the same switch features.

**What does not:** the "concentrates specifically in the GO-invisible, DTU-without-DGE
modules, homogeneously across tissues" refinement. GO-invisible (1.068) and GO-visible
(1.084) are now indistinguishable. **Any drafted prose asserting the GO-invisible
localisation must be rewritten** — this includes `MANUSCRIPT_PLAN.md` §15/§20 and
`GENETIC_ANCHORING_RESULTS.md`, which still carry the old 1.172 / I²=0.00 framing.

### Fig 4A — the SNCA panel is illustrative, not evidence

`06_switch_mechanism/_m/switch_orthogonal_confirm/ORTHOGONAL_CONFIRMATION.md` records that
SNCA's anchored isoform sits at **0.29% of the gene's long-read output**
(`max_anchored_if = 0.0029`, `confirmed_at_usable_abundance = FALSE`) and states plainly
that no single locus should carry a main figure on its own. The panel is retained as a
mechanism vignette and the legend must say so; the set-level evidence is panels B–D plus
S-real-8. Do not present SNCA as the paper's genetic-anchoring evidence.

Rationale for five mains: Fig 1 establishes the method works where truth is known; Fig 2
establishes the real-data modules are reproducible (answering the "fine-grained partition =
noise?" objection); Fig 3 is the genetic anchoring with the matched-baseline null as its
internal control; Fig 4 resolves variants onto switches; Fig 5 answers the composition
confound. The three-baseline rates (Fig S-real-1) deliberately sit in the supplement
because their job is to **bound** the claim (IsoGraph is not globally superior), not to
advance it.

**Display-item budget (Cell Genomics ≤7, figures + tables).** Fig 1–5 + Table 1
(QTL contrast) = 6. Table 2 (splicing-led genes) would make 7 — at the cap. Recommendation:
keep Table 2 in the supplement, where S8/S9 already carry the same content, and hold a
slot in reserve.

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
| S-real-1 | `manuscript/_m/figures/figBaselineRates.{pdf,png}` | **Bounds the claim:** per-module phenotype rate is driven by switch features (both switch-fed methods win: wgcna_switch_only 0.336 > isograph 0.274 > wgcna_multiplex 0.189 ≈ wgcna_gene 0.180) and GO enrichment by abundance — IsoGraph is not globally superior. |
| S-real-2 | `manuscript/_m/figures/figGwasResolution.{pdf,png}` | MAGMA module-GWAS enrichment is size-confounded; at canonical resolution 5.0 the giant-module artifact disappears (0/8 significant modules are giant vs 18/43 at res 2.0 and 79/99 for gene-level WGCNA) yet the schizophrenia signal survives across six regions. Built by `manuscript/_h/gwas_resolution_figure.R`; summary `05_genetic_anchoring/_m/gwas/GWAS_RESOLUTION_SUMMARY.md`. |
| S-real-3 | `manuscript/_m/figures/figGoInvisible.{pdf,png}` | The GO-invisible SCZD switch modules carry functionally-consequential isoform switching comparable to the genome-wide background, and nearly every member carries a real anticorrelated transcript pair — GO-invisibility is GO's gene-level bias, not low module quality. **Note:** this is a claim about module *content*, and it is unaffected by the Fig 3 refresh; it is now the paper's only support for the GO-invisible framing, since the QTL contrast no longer separates GO-invisible from GO-visible. Built by `manuscript/_h/go_invisible_figure.R`. |
| S-real-4 | `manuscript/_m/figures/figSeparation.{pdf,png}` | **Abundance and isoform structure are separable and the separation adds information:** IsoGraph's per-gene abundance and switch axes are largely orthogonal (median \|r\|≈0.13, 41% of genes \|r\|<0.1, **A**); the de-confounded incremental test finds composition-unique genes in most cohort/regions whose switch channel carries phenotype signal total abundance misses (**C**); e.g. NREP (`ENSG00000134986`) has flat total abundance across diagnosis (p=0.93) but a significant isoform switch (p=7e-5, **B**). Read alongside Fig 5, which reports how many of these survive cell-type adjustment. Built by `manuscript/_h/abundance_structure_figure.R`. |
| S-real-5 | `manuscript/_m/figures/figSwitchConsequence.{pdf,png}` | The coding consequence of the switch axis is **productive UTR/CDS remodeling, not decay**: across 9 structural classes only UTR-remodeled (1.27×, 10/10 regions) and CDS-remodeled (1.04×, 10/10) are enriched under a within-gene permutation null, while NMD routing, biotype switch and coding-status loss are depleted; the signal is indistinguishable between GO-invisible and GO-visible modules. Built by `manuscript/_h/switch_consequence_figure.R`; summary `06_switch_mechanism/_m/SWITCH_CONSEQUENCE_SUMMARY.md`. |
| S-real-6 | `manuscript/_m/figures/figRbpRegulon.{pdf,png}` | Switch modules carry recurrent RBP regulons — motifs (PPRC1/KHDRBS1/ELAVL4 8/10 regions; RNASEL/RBMS3/RBM14/ELAVL3/CPEB4/ADAR 7/10; neuronal ELAV/CPEB families) enriched in switched exons across regions and modules (**A**, **C**) — and a subset of the nominated factors show measured **binding capacity** at their regulon genes' switched exons over those genes' own constitutive exons (**B**; 31/38 directionally preferential, 25/38 at q ≤ 0.05, led by PCBP2 0.116 and MATR3 0.105, median gap 0.024). **Framing is deliberately narrow:** ENCODE eCLIP is HepG2/K562, so this is capacity at alternative-exon sequence, NOT neuronal occupancy of these regulons; the neuronal-CLIP validation program is an honest null. Built by `manuscript/_h/rbp_regulon_figure.R`; summary `07_rbp_regulation/_m/rbp/RBP_REGULON_SUMMARY.md`. |
| S-real-7 | `manuscript/_m/figures/figClinicalConsequence.{pdf,png}` | Clinical consequence of the switch layer: (A) switch genes are more LoF-constrained than genome-wide in every region (median LOEUF 0.72 vs 0.94; Fisher p≈1e-93); (B) switched exons carry lower ClinVar P/LP density than constitutive exons (ratio 0.18/0.21, robust to CDS-only scope) — expected alternative-exon biology; (C) the colocalized splicing-led genes are themselves constrained. Built by `manuscript/_h/clinical_consequence_figure.R`; summary `06_switch_mechanism/_m/CLINICAL_CONSEQUENCE_META.md`. |
| S-real-8 | `manuscript/_m/figures/figOrthogonalConfirm.{pdf,png}` | **The switches are not short-read quantification artifacts.** On ONT long-read DLPFC (Aguzzoli-Heberle 2024 NBT; Bambu; n=12 — independent platform, lab, cohort and quantifier) the genetically anchored switch pairs are switch-like at 0.453 vs an abundance-matched null of 0.252 (empirical P = 5e-4, n=53 pairs), rising to 0.600 vs 0.304 (P = 5e-4) when restricted to the 30 pairs whose anchored isoform is usably expressed; 10/12 splicing-led genes confirm. **The two genes that fail the abundance qualification — SNCA (IF 0.003) and CTSH (0.004) — are shown, not dropped**, because they are the two the manuscript is most tempted to feature. Built by `manuscript/_h/orthogonal_confirm_figure.R`. |
| S-real-9 | `manuscript/_m/figures/figIsaConcordance.{pdf,png}` | **The switch signal is not an artifact of IsoGraph.** Under satuRn (the DTU test used by IsoformSwitchAnalyzeR v2), IsoGraph's switch genes carry elevated DTU evidence in **10 of 12 testable analyses** (rank-biserial 0.06–0.44; GTEx hippocampus P = 2e-99, cortex P = 1e-60), with 2 null (n. accumbens, spinal cord). **The remaining 5 of 17 analyses are vacuous by construction** — IsoGraph called no switch genes there — and are drawn as explicit not-testable rows so the denominator cannot be misread. The continuous rank-based statistic is the primary test: the binary empirical-FDR overlap is a conservative lower bound, because Efron's empirical null assumes a sparse alternative and pervasive brain DTU violates it. Built by `manuscript/_h/isa_concordance_figure.R`. |
| S-real-10 | `manuscript/_m/figures/figColocConvergence.{pdf,png}` | **Colocalizing switch genes do NOT concentrate in particular modules — in any trait.** Extends the SCZ-only Fig 4E content to all five traits (AD, PD, LBD, ALS, SCZ) on both aging partitions, and adds the two comparisons a count panel cannot make. (A) Against a size-matched null that holds module sizes fixed, no (trait, source) cell is more concentrated than chance (permutation P = 0.19–1.00, 10/10 cells). (B) Per-module counts, with modules too small in the tested pool to be testable drawn as open symbols; no module survives FDR (min q = 0.23 over 46 testable rows), and every colocalizing gene sits at its own locus (0/58 rows are multi-gene single-locus). (C) **Why the published SCZ convergence does not survive.** `scz_age_projection` reports coloc genes enriched in MAGMA-anchored modules at 15/31 = 48% vs a 25% background (hypergeometric P = 0.004), using ALL module genes as the denominator. But a gene can only colocalize if it was CLPP-tested, and anchored modules are MAGMA-enriched for the same GWAS that decides which genes enter the test — so the tested pool is already 45% anchored, and against it 48% is null (P = 0.40). The same inflation appears for AD (P = 0.043 -> 0.17). This is ascertainment, not convergence. Built by `manuscript/_h/coloc_convergence_figure.R`. |

### Retired / folded display items

- **`figSczConvergence`** — was built and committed but carried no figure number and could
  not be cited. Folded into **Fig 4E**. `manuscript/_h/scz_convergence_figure.R` is kept
  because it is the reproducible standalone builder for the same result; its per-module
  candidate RBP regulators, which do not fit legibly at half width, are in Table S18.

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
candidate regulators (S18, backs Fig 4E). The synthetic benchmark summary —
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
3. **Two orphans resolved.** `figSczConvergence` → Fig 4E; `figS12` numbered.
4. **Fig 1 rebuilt schematic-first** (`figConceptOverview`) — the paper no longer opens
   with a 24-panel boxplot grid, and it now defines its central object before using it.
5. **Three display items built** from analyses that had none: Fig 5
   (`figCompositionRobustness`), S-real-8 (`figOrthogonalConfirm`), S-real-9
   (`figIsaConcordance`), plus the eCLIP panel folded into S-real-6.
6. **Six supplementary tables added** (S13–S18).
7. **Two overclaiming annotations corrected.** Fig 2A now shows both trusted percentages
   (89% vs 88%) rather than implying a difference the counts do not contain; Fig 4A's
   SNCA panel is documented as illustrative against the repository's own long-read result.

**Still open:**

1. **Graphical abstract.** Fig 1A–C is the natural basis for one, but a Cell Press
   graphical abstract is a separate artwork and has not been drawn. Whether Cell Genomics
   research Articles require one is still unverified (`MANUSCRIPT_PLAN.md` §11).
2. **Prose must be rewritten for the Fig 3 refresh.** `MANUSCRIPT_PLAN.md` §15/§20 and
   `GENETIC_ANCHORING_RESULTS.md` still assert the GO-invisible localisation (1.172,
   I²=0.00) that the refreshed inputs do not support.
3. **`CELL_GENOMICS_IMPACT_PROGRESS.md` Item 2 is stale.** It reports 17/39
   binding-supported RBPs at a median gap of 0.013; the committed
   `rbp_binding_support.parquet` now holds 38 RBPs, 25 supported, median gap 0.024
   (the GC-background rescan). S-real-6 and Table S17 use the parquet.
