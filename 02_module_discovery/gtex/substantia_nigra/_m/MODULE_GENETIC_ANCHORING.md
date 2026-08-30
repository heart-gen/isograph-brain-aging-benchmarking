# Module-level genetic anchoring — gtex-aging/substantia_nigra vs GTEx Brain_Substantia_nigra xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region substantia_nigra`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 719 | True | 0.69 | 0.757 | -0.092 | -0.003 | 0.735 |
| M003 | 357 | False | 0.766 | 1.143 | -0.401 | -0.005 | 0.972 |
| M006 | 214 | True | 0.733 | 0.926 | -0.234 | -0.003 | 0.827 |
| M008 | 168 | True | 1.614 | 1.483 | 0.085 | -0.016 | 0.366 |
| M009 | 147 | False | 0.0 | 0.573 | -37.288 | -0.004 | 1.000 |
| M015 | 86 | False | 0.916 | 0.513 | 0.58 | -0.039 | 0.058 |
| M018 | 72 | True | 1.498 | 1.275 | 0.161 | -0.042 | 0.338 |
| M021 | 58 | False | 0.529 | 0.732 | -0.326 | -0.089 | 0.691 |
| M022 | 52 | True | 0.856 | 0.842 | 0.016 | -0.12 | 0.429 |
| M023 | 49 | True | 0.526 | 1.021 | -0.663 | -0.081 | 0.870 |
| M026 | 41 | True | 0.86 | 0.477 | 0.589 | -0.098 | 0.137 |
| M029 | 30 | True | 0.683 | 0.551 | 0.214 | -0.998 | 0.336 |
| M031 | 26 | True | 2.572 | 0.814 | 1.15 | -0.815 | 0.044 |
| M033 | 25 | True | 0.698 | 0.382 | 0.604 | -2.588 | 0.174 |
| M036 | 22 | True | 1.458 | 0.472 | 1.129 | -1.799 | 0.062 |
| M039 | 21 | True |  |  |  |  | NA |
| M041 | 20 | True | 0.637 | 0.815 | -0.246 | -2.809 | 0.578 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._