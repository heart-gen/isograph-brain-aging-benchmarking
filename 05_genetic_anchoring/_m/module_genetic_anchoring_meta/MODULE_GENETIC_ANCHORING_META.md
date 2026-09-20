# Module-level genetic anchoring — cross-analysis meta

Per phenotype-associated module, the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) vs a size-matched permutation null. Pooled over 17 analyses (SCZD + 3 BrainSEQ aging + 13 GTEx aging).

## By stratum

| stratum | n_modules | n_pos_contrast | n_anchored_p05 | median_contrast | frac_pos |
| --- | --- | --- | --- | --- | --- |
| all_phenosig | 54.0 | 29.0 | 5.0 | 0.013 | 0.54 |
| go_invisible | 34.0 | 18.0 | 0.0 | 0.007 | 0.53 |
| go_visible | 20.0 | 11.0 | 5.0 | 0.034 | 0.55 |
| neurodegen_gtex | 45.0 | 25.0 | 4.0 | 0.049 | 0.56 |
| brainseq_sczd | 6.0 | 3.0 | 1.0 | 0.017 | 0.5 |

## Most-anchored modules (smallest perm_p)

| analysis | region | module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gtex-aging | amygdala | M007 | 156 | False | 2.507 | 0.819 | 1.119 | 0.001 |
| gtex-aging | frontal_cortex_ba9 | M015 | 45 | False | 3.184 | 1.029 | 1.13 | 0.003 |
| brainseq-sczd |  | M005 | 382 | False | 1.522 | 0.963 | 0.458 | 0.006 |
| gtex-aging | cerebellum | M007 | 234 | False | 0.937 | 0.617 | 0.418 | 0.018 |
| gtex-aging | anterior_cingulate_cortex_ba24 | M003 | 630 | False | 1.275 | 0.988 | 0.254 | 0.029 |
| gtex-aging | frontal_cortex_ba9 | M002 | 1034 | False | 1.21 | 1.029 | 0.162 | 0.054 |
| gtex-aging | hypothalamus | M024 | 28 | True | 1.954 | 0.875 | 0.803 | 0.063 |
| gtex-aging | amygdala | M030 | 20 | True | 3.01 | 0.994 | 1.108 | 0.065 |
| gtex-aging | anterior_cingulate_cortex_ba24 | M007 | 234 | True | 1.228 | 0.894 | 0.317 | 0.068 |
| gtex-aging | hypothalamus | M025 | 23 | True | 1.75 | 0.71 | 0.902 | 0.07 |
| gtex-aging | hippocampus | M020 | 53 | False | 1.19 | 0.721 | 0.501 | 0.102 |
| brainseq-sczd |  | M004 | 432 | True | 1.111 | 0.943 | 0.164 | 0.148 |
| gtex-aging | frontal_cortex_ba9 | M004 | 760 | False | 1.215 | 1.082 | 0.116 | 0.151 |
| gtex-aging | amygdala | M020 | 37 | True | 1.626 | 0.954 | 0.533 | 0.162 |
| gtex-aging | hippocampus | M001 | 1362 | False | 0.774 | 0.709 | 0.088 | 0.181 |
| gtex-aging | anterior_cingulate_cortex_ba24 | M024 | 21 | True | 1.186 | 0.677 | 0.56 | 0.187 |
| gtex-aging | anterior_cingulate_cortex_ba24 | M001 | 1458 | True | 0.958 | 0.893 | 0.07 | 0.239 |
| brainseq-sczd |  | M009 | 170 | True | 1.348 | 1.17 | 0.142 | 0.268 |
| gtex-aging | hippocampus | M018 | 57 | True | 1.354 | 1.077 | 0.229 | 0.281 |
| gtex-aging | amygdala | M006 | 211 | True | 1.04 | 0.913 | 0.13 | 0.29 |

_A module with positive contrast + small perm_p is anchored to splicing genetics as a unit — evidence the *program*, not just pooled member genes, is genetically anchored. SCZD is expected to be weak (disease is eQTL-led); the splicing-anchored modules should concentrate in GO-invisible / neurodegen-relevant GTEx analyses._