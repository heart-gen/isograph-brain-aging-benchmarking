# Module-level genetic anchoring — cross-analysis meta

Per phenotype-associated module, the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) vs a size-matched permutation null. Pooled over 17 analyses (SCZD + 3 BrainSEQ aging + 13 GTEx aging).

## By stratum

| stratum | n_modules | n_pos_contrast | n_anchored_p05 | median_contrast | frac_pos |
| --- | --- | --- | --- | --- | --- |
| all_phenosig | 97.0 | 51.0 | 7.0 | 0.016 | 0.53 |
| go_invisible | 67.0 | 39.0 | 4.0 | 0.034 | 0.58 |
| go_visible | 30.0 | 12.0 | 3.0 | -0.053 | 0.4 |
| neurodegen_gtex | 83.0 | 45.0 | 6.0 | 0.016 | 0.54 |
| brainseq_sczd | 10.0 | 5.0 | 1.0 | -0.038 | 0.5 |

## Most-anchored modules (smallest perm_p)

| analysis | region | module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gtex-aging | amygdala | M007 | 137 | False | 2.384 | 0.85 | 1.031 | 0.001 |
| gtex-aging | putamen_basal_ganglia | M009 | 256 | False | 0.837 | 0.482 | 0.552 | 0.006 |
| gtex-aging | cerebellum | M006 | 225 | False | 0.869 | 0.579 | 0.406 | 0.025 |
| gtex-aging | anterior_cingulate_cortex_ba24 | M007 | 131 | True | 1.719 | 1.08 | 0.464 | 0.041 |
| gtex-aging | substantia_nigra | M031 | 26 | True | 2.572 | 0.814 | 1.15 | 0.044 |
| brainseq-sczd |  | M030 | 22 | True | 1.935 | 0.651 | 1.089 | 0.046 |
| gtex-aging | putamen_basal_ganglia | M019 | 92 | True | 1.814 | 1.083 | 0.516 | 0.048 |
| gtex-aging | substantia_nigra | M015 | 86 | False | 0.916 | 0.513 | 0.58 | 0.058 |
| brainseq-sczd |  | M017 | 59 | True | 1.287 | 0.693 | 0.62 | 0.062 |
| gtex-aging | substantia_nigra | M036 | 22 | True | 1.458 | 0.472 | 1.129 | 0.062 |
| gtex-aging | frontal_cortex_ba9 | M002 | 604 | False | 1.083 | 0.913 | 0.172 | 0.082 |
| gtex-aging | hippocampus | M000 | 1369 | False | 0.826 | 0.713 | 0.148 | 0.091 |
| gtex-aging | amygdala | M016 | 31 | True | 1.966 | 0.907 | 0.774 | 0.091 |
| gtex-aging | amygdala | M008 | 115 | True | 1.147 | 0.741 | 0.437 | 0.096 |
| gtex-aging | frontal_cortex_ba9 | M022 | 23 | True | 2.423 | 1.149 | 0.746 | 0.102 |
| gtex-aging | frontal_cortex_ba9 | M014 | 81 | True | 1.287 | 0.862 | 0.401 | 0.112 |
| gtex-aging | amygdala | M018 | 26 | True | 0.373 | 0.17 | 0.787 | 0.115 |
| gtex-aging | substantia_nigra | M026 | 41 | True | 0.86 | 0.477 | 0.589 | 0.137 |
| gtex-aging | hippocampus | M004 | 272 | False | 1.346 | 1.102 | 0.201 | 0.148 |
| brainseq-aging | caudate | M008 | 205 | True | 1.105 | 0.881 | 0.226 | 0.172 |

_A module with positive contrast + small perm_p is anchored to splicing genetics as a unit — evidence the *program*, not just pooled member genes, is genetically anchored. SCZD is expected to be weak (disease is eQTL-led); the splicing-anchored modules should concentrate in GO-invisible / neurodegen-relevant GTEx analyses._