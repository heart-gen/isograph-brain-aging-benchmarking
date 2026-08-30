# Module-level genetic anchoring — gtex-aging/anterior_cingulate_cortex_ba24 vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 550 | False | 0.488 | 0.72 | -0.389 | -0.007 | 0.994 |
| M001 | 513 | False | 0.96 | 0.866 | 0.104 | 0.001 | 0.253 |
| M002 | 456 | False | 1.024 | 1.165 | -0.129 | -0.001 | 0.781 |
| M004 | 320 | True | 0.756 | 0.794 | -0.048 | 0.004 | 0.633 |
| M006 | 280 | True | 1.103 | 1.021 | 0.076 | -0.006 | 0.349 |
| M007 | 263 | True | 1.529 | 1.556 | -0.018 | 0.0 | 0.540 |
| M008 | 255 | False | 0.914 | 0.776 | 0.163 | -0.01 | 0.214 |
| M010 | 178 | True | 0.827 | 1.091 | -0.278 | -0.013 | 0.856 |
| M011 | 165 | False | 2.151 | 0.71 | 1.109 | -0.029 | 0.001 |
| M013 | 155 | True | 1.121 | 0.891 | 0.229 | -0.024 | 0.176 |
| M015 | 142 | True | 1.461 | 1.141 | 0.248 | -0.027 | 0.167 |
| M016 | 121 | True | 0.824 | 1.419 | -0.543 | -0.028 | 0.945 |
| M017 | 115 | True | 1.213 | 0.997 | 0.196 | -0.018 | 0.246 |
| M019 | 100 | True | 0.559 | 0.755 | -0.301 | -0.018 | 0.793 |
| M020 | 85 | True | 0.718 | 0.852 | -0.17 | -0.039 | 0.653 |
| M021 | 80 | True | 0.582 | 1.477 | -0.931 | -0.024 | 0.987 |
| M023 | 77 | False | 0.721 | 0.586 | 0.207 | -0.041 | 0.269 |
| M027 | 52 | True | 2.338 | 1.683 | 0.329 | -0.018 | 0.253 |
| M028 | 43 | True | 0.944 | 1.476 | -0.447 | -0.078 | 0.765 |
| M031 | 29 | False | 0.637 | 1.12 | -0.564 | -0.326 | 0.791 |
| M032 | 29 | True | 2.864 | 0.761 | 1.325 | -0.163 | 0.011 |
| M033 | 29 | True | 0.709 | 0.983 | -0.326 | -0.797 | 0.658 |
| M035 | 26 | True | 1.186 | 1.292 | -0.085 | -0.195 | 0.514 |
| M039 | 21 | True | 0.888 | 1.822 | -0.719 | -1.277 | 0.792 |
| M040 | 21 | True | 0.809 | 0.438 | 0.614 | -0.852 | 0.188 |
| M041 | 21 | True | 1.574 | 1.8 | -0.134 | -0.925 | 0.536 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._