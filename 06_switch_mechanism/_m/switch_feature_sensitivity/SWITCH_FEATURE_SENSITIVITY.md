# Preprocessing sensitivity of the switch representation — all regions

`switch_feature_sensitivity.py --aggregate` over **14 cohort × region fits**. Per-region tables are in `<cohort>/<region>/SWITCH_FEATURE_SENSITIVITY.md`.

## Scope

The module partition is held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. It measures the stability of the representation and its trait signal, not of an independently refit network. Every region passed the exact-rebuild gate (largest max |diff| = 1.18e-13). Published settings are cohort-specific — BrainSEQ is fit on transcripts with count > 10 in ≥ 70% of samples, GTEx on the unfiltered matrix — so the GTEx expression axis runs from no filter to the strictest filter.

## 1-3. Worst case over every non-published setting, per region

`min_retained_frac` is the smallest fraction of the region's published FDR-significant module–age associations that survive any single setting (blank where the region has none); `max_sign_flips_among_published_sig` counts reversals of a reported effect.

| cohort | region | n_settings | n_fdr_sig_published | min_effect_pearson | min_retained_frac | max_sign_flips_among_published_sig | min_median_abs_feature_r |
|---|---|---|---|---|---|---|---|
| brainseq | caudate | 11 | 20 | 0.889 | 0.75 | 0 | 0.842 |
| gtex | amygdala | 12 | 11 | 0.207 | 0.0909 | 3 | 0.228 |
| gtex | anterior_cingulate_cortex_ba24 | 12 | 20 | 0.687 | 0.65 | 1 | 0.253 |
| gtex | caudate_basal_ganglia | 12 | 0 | 0.404 |  | 0 | 0.23 |
| gtex | cerebellar_hemisphere | 12 | 12 | 0.295 | 0.5 | 5 | 0.286 |
| gtex | cerebellum | 12 | 7 | 0.429 | 0.857 | 1 | 0.19 |
| gtex | cortex | 12 | 24 | 0.341 | 0.25 | 8 | 0.197 |
| gtex | frontal_cortex_ba9 | 12 | 40 | 0.583 | 0.55 | 12 | 0.252 |
| gtex | hippocampus | 12 | 17 | 0.534 | 0.588 | 4 | 0.236 |
| gtex | hypothalamus | 12 | 8 | 0.119 | 0.125 | 5 | 0.248 |
| gtex | nucleus_accumbens_basal_ganglia | 12 | 0 | 0.0702 |  | 0 | 0.291 |
| gtex | putamen_basal_ganglia | 12 | 0 | 0.502 |  | 0 | 0.235 |
| gtex | spinal_cord_cervical_c_1 | 12 | 0 | 0.468 |  | 0 | 0.219 |
| gtex | substantia_nigra | 12 | 9 | 0.289 | 0.556 | 1 | 0.225 |

Minimum module age-effect correlation with the published effects, by axis:

| cohort | region | expression | minor_isoform | pseudocount |
|---|---|---|---|---|
| brainseq | caudate | 0.889 | 0.896 | 0.992 |
| gtex | amygdala | 0.207 | 0.375 | 0.856 |
| gtex | anterior_cingulate_cortex_ba24 | 0.687 | 0.804 | 0.955 |
| gtex | caudate_basal_ganglia | 0.404 | 0.551 | 0.825 |
| gtex | cerebellar_hemisphere | 0.295 | 0.503 | 0.982 |
| gtex | cerebellum | 0.429 | 0.467 | 0.975 |
| gtex | cortex | 0.341 | 0.504 | 0.975 |
| gtex | frontal_cortex_ba9 | 0.583 | 0.731 | 0.976 |
| gtex | hippocampus | 0.534 | 0.554 | 0.832 |
| gtex | hypothalamus | 0.119 | 0.22 | 0.81 |
| gtex | nucleus_accumbens_basal_ganglia | 0.0702 | 0.252 | 0.84 |
| gtex | putamen_basal_ganglia | 0.502 | 0.558 | 0.821 |
| gtex | spinal_cord_cervical_c_1 | 0.468 | 0.602 | 0.847 |
| gtex | substantia_nigra | 0.289 | 0.487 | 0.87 |

## 4. Identifiability — median |switch–age r| by transcript-number stratum

| cohort | region | [10, 20) | [2, 3) | [20, 10000) | [3, 5) | [5, 10) |
|---|---|---|---|---|---|---|
| brainseq | caudate | 0.0842 | 0.0657 | 0.102 | 0.0691 | 0.0714 |
| gtex | amygdala | 0.0731 | 0.073 | 0.079 | 0.0665 | 0.066 |
| gtex | anterior_cingulate_cortex_ba24 | 0.0732 | 0.0755 | 0.0873 | 0.0661 | 0.0649 |
| gtex | caudate_basal_ganglia | 0.0479 | 0.0499 | 0.0481 | 0.0452 | 0.0466 |
| gtex | cerebellar_hemisphere | 0.0784 | 0.0644 | 0.102 | 0.0642 | 0.0662 |
| gtex | cerebellum | 0.0573 | 0.0578 | 0.0645 | 0.0563 | 0.0535 |
| gtex | cortex | 0.0672 | 0.0714 | 0.0801 | 0.0625 | 0.0608 |
| gtex | frontal_cortex_ba9 | 0.0799 | 0.0726 | 0.0998 | 0.0661 | 0.0704 |
| gtex | hippocampus | 0.0726 | 0.0703 | 0.0872 | 0.0617 | 0.0666 |
| gtex | hypothalamus | 0.0686 | 0.0639 | 0.0794 | 0.0597 | 0.0646 |
| gtex | nucleus_accumbens_basal_ganglia | 0.0669 | 0.0578 | 0.077 | 0.0578 | 0.0572 |
| gtex | putamen_basal_ganglia | 0.0564 | 0.0548 | 0.0638 | 0.0531 | 0.0542 |
| gtex | spinal_cord_cervical_c_1 | 0.0476 | 0.0526 | 0.0477 | 0.0492 | 0.0485 |
| gtex | substantia_nigra | 0.0763 | 0.0734 | 0.0889 | 0.0678 | 0.0696 |

## 5. Quantification pipeline, read against a same-quantifier reference

`cross_quantifier` rows compare matched regions across cohorts (Salmon vs RSEM), which confounds quantifier with cohort. `same_quantifier` rows hold the quantifier fixed and vary the tissue within a cohort, usually on overlapping donors (`n_shared_donors`), which inflates their concordance through shared donor effects. A cross-cohort figure is only evidence about the quantifier to the extent it falls below the same-quantifier reference.

| pair_type | cohort_1 | region_1 | quantifier_1 | cohort_2 | region_2 | quantifier_2 | n_samples_1 | n_samples_2 | n_shared_donors | n_shared_genes | pearson_gene_age_effect | spearman_gene_age_effect | sign_concordance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cross_quantifier | brainseq | caudate | Salmon | gtex | caudate_basal_ganglia | RSEM | 238 | 300 | 0 | 11476 | 0.00707 | 0.00309 | 0.499 |
| cross_quantifier | brainseq | hippocampus | Salmon | gtex | hippocampus | RSEM | 238 | 255 | 0 | 10944 | -0.00159 | 0.00206 | 0.499 |
| cross_quantifier | brainseq | dlpfc | Salmon | gtex | frontal_cortex_ba9 | RSEM | 222 | 269 | 0 | 11340 | 0.000357 | 0.00157 | 0.509 |
| same_quantifier | gtex | caudate_basal_ganglia | RSEM | gtex | putamen_basal_ganglia | RSEM | 300 | 254 | 233 | 16073 | 0.326 | 0.308 | 0.606 |
| same_quantifier | gtex | caudate_basal_ganglia | RSEM | gtex | nucleus_accumbens_basal_ganglia | RSEM | 300 | 285 | 249 | 16203 | 0.306 | 0.296 | 0.6 |
| same_quantifier | gtex | frontal_cortex_ba9 | RSEM | gtex | cortex | RSEM | 269 | 270 | 194 | 16222 | 0.514 | 0.474 | 0.658 |
| same_quantifier | gtex | cerebellum | RSEM | gtex | cerebellar_hemisphere | RSEM | 266 | 277 | 222 | 16062 | 0.402 | 0.384 | 0.633 |
| same_quantifier | brainseq | caudate | Salmon | brainseq | dlpfc | Salmon | 238 | 222 | 173 | 10840 | 0.227 | 0.206 | 0.566 |
| same_quantifier | brainseq | caudate | Salmon | brainseq | hippocampus | Salmon | 238 | 238 | 199 | 10629 | 0.175 | 0.176 | 0.559 |
| same_quantifier | brainseq | dlpfc | Salmon | brainseq | hippocampus | Salmon | 222 | 238 | 188 | 10845 | 0.201 | 0.18 | 0.557 |
