# Manuscript figure & table ordering (IsoGraph, Nature Methods)

The final figure/table sequence for the IsoGraph manuscript, with the one honest claim each
artifact carries. Synthetic-benchmark figures live under `benchmark/03_metrics/figures/`;
real-data figures under `real_data/{_m,stability/_m,gwas/_m}/figures/`. North-star: IsoGraph is
a **complementary DTU-without-DGE layer**, not a globally superior method — the ordering moves
from "the method recovers switch modules" (synthetic) → "its modules are trustworthy"
(real-data reproducibility) → "they carry real, genetically-anchored biology invisible to
abundance pipelines" (real-data complementarity).

## Main figures

| Fig | File | Claim |
| --- | --- | --- |
| **1** | `benchmark/03_metrics/figures/fig1_benchmark_overview.{pdf,png}` | On synthetic ground truth IsoGraph VAE recovers switch modules with complete switch-gene detection; 219/240 paired Wilcoxon tests favour it over WGCNA. |
| **2** | `real_data/stability/_m/figures/figTrustFunnel.{pdf,png}` | On real brain data IsoGraph's modules are per-module trustworthy: 236/266 chance-trusted across six regions, switch drivers reproduce (ρ≈0.77–0.82), 25 modules replicate aging cross-cohort (~4× the abundance baseline). |
| **3** | `real_data/_m/figures/figQtlSpecificity.{pdf,png}` | The disease/GO-invisible switch modules are genetically anchored — splicing QTL are spared relative to eQTL (ratio≈1.13) exactly there and nowhere in the GO-visible control, and the effect is IsoGraph-only on matched WGCNA baselines. |

Rationale for three mains: Fig 1 establishes the method works where truth is known; Fig 2
establishes the real-data modules are reproducible (answering the "fine-grained partition =
noise?" objection); Fig 3 is the payoff — orthogonal genetic evidence that the complementary
layer is real, with two internal controls (GO-visible null, matched-baseline null). The
three-baseline rates (Fig S-real-1) deliberately sit in the supplement because their job is to
**bound** the claim (IsoGraph is not globally superior), not to advance it.

## Supplementary figures

Synthetic parameter-resolved and stress panels (existing, `benchmark/03_metrics/figures/`):

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

Real-data supplements:

| Fig | File | Claim |
| --- | --- | --- |
| S-real-1 | `real_data/_m/figures/figBaselineRates.{pdf,png}` | **Bounds the claim:** per-module phenotype rate is driven by switch features (both switch-fed methods win) and GO enrichment by abundance — IsoGraph is not globally superior. |
| S-real-2 | `real_data/_m/figures/figGwasResolution.{pdf,png}` | MAGMA module-GWAS enrichment is size-confounded; at canonical resolution 5.0 the giant-module artifact disappears (0/8 significant modules are giant vs 18/43 at res 2.0 and 79/99 for gene-level WGCNA) yet the schizophrenia signal survives across six regions. Built to manuscript conventions by `real_data/_h/gwas_resolution_figure.R`; summary `real_data/gwas/_m/GWAS_RESOLUTION_SUMMARY.md`. |
| S-real-3 | `real_data/_m/figures/figGoInvisible.{pdf,png}` | The four GO-invisible SCZD switch modules carry functionally-consequential isoform switching comparable to the genome-wide background, and nearly every member carries a real anticorrelated transcript pair — GO-invisibility is GO's gene-level bias, not low module quality. Built by `real_data/_h/go_invisible_figure.R`; sits beside Fig 3. |
| S-real-4 | `real_data/_m/figures/figSeparation.{pdf,png}` | **Abundance and isoform structure are separable and the separation adds information:** IsoGraph's per-gene abundance and switch axes are largely orthogonal (median \|r\|≈0.13, 41% of genes \|r\|<0.1, **A**); the de-confounded incremental test finds a specific set of composition-unique genes in every cohort/region whose switch channel carries phenotype signal total abundance misses (34 SCZD / 43 aging-caudate, up to 545 in GTEx cortex; **C**); e.g. NREP (`ENSG00000134986`) has flat total abundance across diagnosis (p=0.93) but a significant isoform switch (p=7e-5, **B**). Built by `real_data/_h/abundance_structure_figure.R` from `abundance_structure_separation.py` outputs; supports the complementary-layer framing (this is non-redundancy, not superiority). |

## Supplementary tables

See `real_data/_m/supp_tables/SUPPLEMENTARY_TABLES.md` (Tables S1–S7) — three-baseline pooled
(S1) and per-region (S2); QTL specificity contrast (S3), matched-baseline contrast (S4) and
raw cis-QTL ORs (S5); GO-invisible disease modules (S6); per-region trust funnel (S7). The
synthetic Table 1 (`benchmark/03_metrics/_m/table1_benchmark_summary.csv`) and scale-compute
table (`tableS_scale_compute_summary.csv`) accompany Figs 1/S6.

## Status

Both previously-open gaps are now closed:

1. **GWAS res-5.0 figure (S-real-2)** — rebuilt to manuscript-figure conventions
   (`figGwasResolution`, Okabe-Ito, no in-panel titles) and given a Manubot summary
   (`real_data/gwas/_m/GWAS_RESOLUTION_SUMMARY.md`). The legacy MAGMA dotplots
   (`real_data/gwas/_m/figures/magma_enrichment_*`) are superseded for the manuscript.
2. **GO-invisible functional-consequence panel (S-real-3)** — built (`figGoInvisible`) to sit
   beside Fig 3, showing the four disease modules' consequence fractions against background
   and their switch coherence.

All main (Fig 1–3) and supplementary (S1–S11, S-real-1/2/3/4) figures now exist; the real-data
figures were generated this session and are not yet committed. S-real-4 (`figSeparation`) states
the foundational non-redundancy claim directly — abundance and isoform structure are
computationally separable and the separation adds phenotype information.
