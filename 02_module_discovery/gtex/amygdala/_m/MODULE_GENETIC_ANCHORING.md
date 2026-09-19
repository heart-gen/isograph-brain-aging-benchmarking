# Module-level genetic anchoring — gtex-aging/amygdala vs GTEx Brain_Amygdala xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region amygdala`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 1370 | False | 0.799 | 0.945 | -0.168 | 0.005 | 0.945 |
| M001 | 1252 | False | 0.647 | 0.76 | -0.161 | -0.006 | 0.927 |
| M006 | 211 | True | 1.04 | 0.913 | 0.13 | -0.01 | 0.290 |
| M007 | 156 | False | 2.507 | 0.819 | 1.119 | -0.013 | 0.001 |
| M010 | 88 | True | 0.995 | 0.906 | 0.094 | -0.037 | 0.392 |
| M012 | 68 | False | 0.448 | 0.703 | -0.451 | -0.061 | 0.830 |
| M020 | 37 | True | 1.626 | 0.954 | 0.533 | -0.127 | 0.162 |
| M022 | 30 | True | 1.433 | 1.365 | 0.049 | -1.176 | 0.456 |
| M030 | 20 | True | 3.01 | 0.994 | 1.108 | -3.342 | 0.065 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._