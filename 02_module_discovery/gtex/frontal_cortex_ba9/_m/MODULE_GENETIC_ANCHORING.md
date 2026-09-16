# Module-level genetic anchoring — gtex-aging/frontal_cortex_ba9 vs GTEx Brain_Frontal_Cortex_BA9 xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 2467 | False | 0.937 | 1.089 | -0.15 | -0.001 | 0.968 |
| M001 | 1187 | True | 0.74 | 0.903 | -0.198 | 0.0 | 0.975 |
| M002 | 604 | False | 1.083 | 0.913 | 0.172 | -0.006 | 0.082 |
| M004 | 448 | True | 1.24 | 1.095 | 0.124 | -0.006 | 0.204 |
| M005 | 178 | True | 0.961 | 0.929 | 0.034 | -0.002 | 0.452 |
| M010 | 92 | True | 1.261 | 1.039 | 0.193 | -0.014 | 0.270 |
| M011 | 92 | True | 1.682 | 1.38 | 0.198 | -0.025 | 0.252 |
| M014 | 81 | True | 1.287 | 0.862 | 0.401 | -0.027 | 0.112 |
| M015 | 79 | True | 1.027 | 1.231 | -0.181 | -0.025 | 0.677 |
| M018 | 59 | True | 1.098 | 0.941 | 0.154 | -0.052 | 0.336 |
| M019 | 35 | True | 1.319 | 1.497 | -0.126 | -0.135 | 0.567 |
| M022 | 23 | True | 2.423 | 1.149 | 0.746 | -0.139 | 0.102 |
| M023 | 22 | True | 1.352 | 1.261 | 0.07 | -0.282 | 0.418 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._