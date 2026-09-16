# Preprocessing sensitivity of the switch representation — all regions

`switch_feature_sensitivity.py --aggregate` over **14 cohort × region fits**. Per-region tables are in `<cohort>/<region>/SWITCH_FEATURE_SENSITIVITY.md`.

## Scope

The module partition is held fixed at the published one; preprocessing is varied and the features, eigengenes and age association are recomputed. It measures the stability of the representation and its trait signal, not of an independently refit network. Every region passed the exact-rebuild gate (largest max |diff| = 0). Both cohorts are fit on the switching transcript filter, so the transcript-filter (`expression`) axis is the same grid in every region.

## 1-3. Worst case over every non-published setting, per region

`min_retained_frac` is the smallest fraction of the region's published FDR-significant module–age associations that survive any single setting (blank where the region has none); `max_sign_flips_among_published_sig` counts reversals of a reported effect.

| cohort | region | n_settings | n_fdr_sig_published | min_effect_pearson | min_retained_frac | max_sign_flips_among_published_sig | min_median_abs_feature_r |
|---|---|---|---|---|---|---|---|
| brainseq | caudate | 13 | 4 | 0.94 | 0.75 | 0 | 0.736 |
| gtex | amygdala | 13 | 6 | 0.781 | 0.667 | 0 | 0.822 |
| gtex | anterior_cingulate_cortex_ba24 | 13 | 4 | 0.64 | 0 | 1 | 0.825 |
| gtex | caudate_basal_ganglia | 13 | 0 | 0.622 |  | 0 | 0.79 |
| gtex | cerebellar_hemisphere | 13 | 6 | 0.765 | 0.833 | 0 | 0.769 |
| gtex | cerebellum | 13 | 11 | 0.726 | 0.818 | 0 | 0.642 |
| gtex | cortex | 13 | 7 | 0.827 | 0.286 | 1 | 0.755 |
| gtex | frontal_cortex_ba9 | 13 | 10 | 0.911 | 0.3 | 1 | 0.811 |
| gtex | hippocampus | 13 | 9 | 0.717 | 0.667 | 2 | 0.817 |
| gtex | hypothalamus | 13 | 8 | 0.731 | 0 | 1 | 0.793 |
| gtex | nucleus_accumbens_basal_ganglia | 13 | 5 | 0.852 | 0 | 1 | 0.821 |
| gtex | putamen_basal_ganglia | 13 | 0 | 0.702 |  | 0 | 0.803 |
| gtex | spinal_cord_cervical_c_1 | 13 | 0 | 0.508 |  | 0 | 0.794 |
| gtex | substantia_nigra | 13 | 0 | 0.621 |  | 0 | 0.79 |

Minimum module age-effect correlation with the published effects, by axis:

| cohort | region | expression | minor_isoform | pseudocount |
|---|---|---|---|---|
| brainseq | caudate | 0.94 | 0.972 | 0.996 |
| gtex | amygdala | 0.781 | 0.82 | 0.981 |
| gtex | anterior_cingulate_cortex_ba24 | 0.64 | 0.909 | 0.979 |
| gtex | caudate_basal_ganglia | 0.622 | 0.879 | 0.973 |
| gtex | cerebellar_hemisphere | 0.765 | 0.87 | 0.992 |
| gtex | cerebellum | 0.726 | 0.849 | 0.993 |
| gtex | cortex | 0.827 | 0.842 | 0.998 |
| gtex | frontal_cortex_ba9 | 0.911 | 0.915 | 0.994 |
| gtex | hippocampus | 0.717 | 0.889 | 0.985 |
| gtex | hypothalamus | 0.731 | 0.831 | 0.987 |
| gtex | nucleus_accumbens_basal_ganglia | 0.852 | 0.888 | 0.99 |
| gtex | putamen_basal_ganglia | 0.702 | 0.797 | 0.992 |
| gtex | spinal_cord_cervical_c_1 | 0.508 | 0.666 | 0.945 |
| gtex | substantia_nigra | 0.621 | 0.753 | 0.983 |

## 4. Identifiability — median |switch–age r| by transcript-number stratum

| cohort | region | [10, 20) | [2, 3) | [3, 5) | [5, 10) |
|---|---|---|---|---|---|
| brainseq | caudate | 0.0654 | 0.0616 | 0.0608 | 0.0568 |
| gtex | amygdala | 0.061 | 0.0788 | 0.0749 | 0.0681 |
| gtex | anterior_cingulate_cortex_ba24 | 0.0776 | 0.0734 | 0.0746 | 0.068 |
| gtex | caudate_basal_ganglia | 0.0658 | 0.0518 | 0.0484 | 0.0454 |
| gtex | cerebellar_hemisphere | 0.0594 | 0.0779 | 0.0776 | 0.0739 |
| gtex | cerebellum | 0.0571 | 0.0673 | 0.0639 | 0.0601 |
| gtex | cortex | 0.0553 | 0.0796 | 0.0737 | 0.0686 |
| gtex | frontal_cortex_ba9 | 0.0606 | 0.0844 | 0.0798 | 0.0724 |
| gtex | hippocampus | 0.0711 | 0.0704 | 0.0675 | 0.0698 |
| gtex | hypothalamus | 0.0401 | 0.0785 | 0.0703 | 0.0677 |
| gtex | nucleus_accumbens_basal_ganglia | 0.0479 | 0.0705 | 0.068 | 0.059 |
| gtex | putamen_basal_ganglia | 0.0479 | 0.0623 | 0.0572 | 0.0551 |
| gtex | spinal_cord_cervical_c_1 | 0.0607 | 0.0442 | 0.0456 | 0.0474 |
| gtex | substantia_nigra | 0.0756 | 0.0722 | 0.0685 | 0.0674 |

## 5. Quantification pipeline, read against a same-quantifier reference

`cross_quantifier` rows compare matched regions across cohorts (Salmon vs RSEM), which confounds quantifier with cohort. `same_quantifier` rows hold the quantifier fixed and vary the tissue within a cohort, usually on overlapping donors (`n_shared_donors`), which inflates their concordance through shared donor effects. A cross-cohort figure is only evidence about the quantifier to the extent it falls below the same-quantifier reference.

| pair_type | cohort_1 | region_1 | quantifier_1 | cohort_2 | region_2 | quantifier_2 | n_samples_1 | n_samples_2 | n_shared_donors | n_shared_genes | pearson_gene_age_effect | spearman_gene_age_effect | sign_concordance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cross_quantifier | brainseq | caudate | Salmon | gtex | caudate_basal_ganglia | RSEM | 238 | 300 | 0 | 10988 | 0.0218 | 0.0184 | 0.511 |
| cross_quantifier | brainseq | hippocampus | Salmon | gtex | hippocampus | RSEM | 238 | 255 | 0 | 10861 | 0.00707 | 0.00586 | 0.506 |
| cross_quantifier | brainseq | dlpfc | Salmon | gtex | frontal_cortex_ba9 | RSEM | 222 | 269 | 0 | 10828 | 0.0184 | 0.0261 | 0.511 |
| same_quantifier | gtex | caudate_basal_ganglia | RSEM | gtex | putamen_basal_ganglia | RSEM | 300 | 254 | 233 | 11946 | 0.348 | 0.335 | 0.616 |
| same_quantifier | gtex | caudate_basal_ganglia | RSEM | gtex | nucleus_accumbens_basal_ganglia | RSEM | 300 | 285 | 249 | 12087 | 0.321 | 0.312 | 0.61 |
| same_quantifier | gtex | frontal_cortex_ba9 | RSEM | gtex | cortex | RSEM | 269 | 270 | 194 | 12037 | 0.588 | 0.553 | 0.688 |
| same_quantifier | gtex | cerebellum | RSEM | gtex | cerebellar_hemisphere | RSEM | 266 | 277 | 222 | 12119 | 0.511 | 0.498 | 0.675 |
| same_quantifier | brainseq | caudate | Salmon | brainseq | dlpfc | Salmon | 238 | 222 | 173 | 12045 | 0.268 | 0.229 | 0.574 |
| same_quantifier | brainseq | caudate | Salmon | brainseq | hippocampus | Salmon | 238 | 238 | 199 | 12207 | 0.207 | 0.197 | 0.567 |
| same_quantifier | brainseq | dlpfc | Salmon | brainseq | hippocampus | Salmon | 222 | 238 | 188 | 12328 | 0.198 | 0.172 | 0.557 |
