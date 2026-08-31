# Module-level genetic anchoring — gtex-aging/amygdala vs GTEx Brain_Amygdala xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region amygdala`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 710 | False | 0.806 | 0.858 | -0.063 | -0.004 | 0.660 |
| M001 | 418 | True | 0.641 | 0.883 | -0.322 | -0.0 | 0.950 |
| M002 | 379 | True | 1.018 | 1.03 | -0.011 | -0.002 | 0.527 |
| M004 | 279 | True | 1.028 | 1.275 | -0.216 | -0.008 | 0.807 |
| M006 | 223 | False | 1.024 | 0.784 | 0.267 | -0.013 | 0.130 |
| M009 | 127 | True | 2.081 | 0.927 | 0.808 | -0.024 | 0.004 |
| M012 | 81 | True | 0.916 | 1.192 | -0.264 | -0.044 | 0.705 |
| M014 | 58 | True | 1.005 | 0.93 | 0.078 | -0.131 | 0.429 |
| M015 | 52 | True | 1.09 | 0.635 | 0.541 | -0.31 | 0.147 |
| M017 | 40 | True | 0.867 | 1.135 | -0.269 | -0.62 | 0.647 |
| M020 | 31 | True | 0.762 | 0.731 | 0.042 | -0.701 | 0.471 |
| M021 | 30 | True | 3.922 | 2.281 | 0.542 | -2.749 | 0.207 |
| M022 | 29 | True | 0.812 | 0.293 | 1.019 | -1.457 | 0.055 |
| M025 | 20 | False | 1.335 | 0.458 | 1.07 | -3.248 | 0.101 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._