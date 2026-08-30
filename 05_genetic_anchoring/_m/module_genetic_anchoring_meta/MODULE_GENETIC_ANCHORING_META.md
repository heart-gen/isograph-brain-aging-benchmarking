# Module-level genetic anchoring — cross-analysis meta

Per phenotype-associated module, the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) vs a size-matched permutation null. Pooled over 17 analyses (SCZD + 3 BrainSEQ aging + 13 GTEx aging).

## By stratum

| stratum | n_modules | n_pos_contrast | n_anchored_p05 | median_contrast | frac_pos |
| --- | --- | --- | --- | --- | --- |
| all_phenosig | 177.0 | 84.0 | 14.0 | -0.034 | 0.47 |
| go_invisible | 124.0 | 57.0 | 7.0 | -0.041 | 0.46 |
| go_visible | 53.0 | 27.0 | 7.0 | 0.008 | 0.51 |
| neurodegen_gtex | 156.0 | 73.0 | 11.0 | -0.034 | 0.47 |
| brainseq_sczd | 7.0 | 4.0 | 1.0 | 0.098 | 0.57 |

## Most-anchored modules (smallest perm_p)

| analysis | region | module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gtex-aging | anterior_cingulate_cortex_ba24 | M011 | 165 | False | 2.151 | 0.71 | 1.109 | 0.001 |
| gtex-aging | hippocampus | M026 | 43 | False | 1.461 | 0.332 | 1.483 | 0.001 |
| brainseq-aging | caudate | M005 | 331 | True | 1.637 | 0.956 | 0.538 | 0.002 |
| gtex-aging | amygdala | M009 | 127 | True | 2.081 | 0.927 | 0.808 | 0.004 |
| gtex-aging | hypothalamus | M005 | 240 | True | 2.127 | 1.257 | 0.526 | 0.004 |
| gtex-aging | putamen_basal_ganglia | M009 | 256 | False | 0.837 | 0.482 | 0.552 | 0.006 |
| gtex-aging | anterior_cingulate_cortex_ba24 | M032 | 29 | True | 2.864 | 0.761 | 1.325 | 0.011 |
| brainseq-sczd |  | M012 | 58 | False | 0.821 | 0.314 | 0.96 | 0.019 |
| gtex-aging | frontal_cortex_ba9 | M005 | 504 | False | 1.676 | 1.25 | 0.293 | 0.02 |
| gtex-aging | cerebellum | M002 | 234 | False | 0.846 | 0.559 | 0.413 | 0.027 |
| brainseq-aging | caudate | M002 | 561 | True | 1.41 | 1.092 | 0.255 | 0.036 |
| gtex-aging | substantia_nigra | M031 | 26 | True | 2.572 | 0.814 | 1.15 | 0.044 |
| gtex-aging | hippocampus | M029 | 35 | False | 2.105 | 0.89 | 0.861 | 0.046 |
| gtex-aging | putamen_basal_ganglia | M019 | 92 | True | 1.814 | 1.083 | 0.516 | 0.048 |
| gtex-aging | hippocampus | M001 | 779 | False | 1.117 | 0.912 | 0.203 | 0.054 |
| gtex-aging | amygdala | M022 | 29 | True | 0.812 | 0.293 | 1.019 | 0.055 |
| gtex-aging | substantia_nigra | M015 | 86 | False | 0.916 | 0.513 | 0.58 | 0.058 |
| gtex-aging | cortex | M005 | 157 | True | 1.184 | 0.811 | 0.378 | 0.061 |
| gtex-aging | substantia_nigra | M036 | 22 | True | 1.458 | 0.472 | 1.129 | 0.062 |
| brainseq-aging | caudate | M042 | 23 | True | 1.756 | 0.646 | 1.0 | 0.064 |

_A module with positive contrast + small perm_p is anchored to splicing genetics as a unit — evidence the *program*, not just pooled member genes, is genetically anchored. SCZD is expected to be weak (disease is eQTL-led); the splicing-anchored modules should concentrate in GO-invisible / neurodegen-relevant GTEx analyses._