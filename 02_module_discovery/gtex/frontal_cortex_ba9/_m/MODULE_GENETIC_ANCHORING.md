# Module-level genetic anchoring — gtex-aging/frontal_cortex_ba9 vs GTEx Brain_Frontal_Cortex_BA9 xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 3685 | False | 0.76 | 0.969 | -0.243 | 0.001 | 0.999 |
| M002 | 1034 | False | 1.21 | 1.029 | 0.162 | 0.003 | 0.054 |
| M004 | 760 | False | 1.215 | 1.082 | 0.116 | -0.006 | 0.151 |
| M005 | 526 | True | 0.99 | 0.981 | 0.01 | -0.001 | 0.462 |
| M008 | 156 | True | 1.041 | 1.062 | -0.02 | -0.005 | 0.538 |
| M010 | 95 | True | 1.058 | 1.187 | -0.115 | -0.05 | 0.586 |
| M015 | 45 | False | 3.184 | 1.029 | 1.13 | -0.061 | 0.003 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._