# Module-level genetic anchoring — gtex-aging/anterior_cingulate_cortex_ba24 vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 2517 | False | 0.79 | 0.922 | -0.155 | -0.001 | 0.964 |
| M001 | 1096 | True | 0.948 | 0.963 | -0.016 | -0.002 | 0.526 |
| M002 | 370 | False | 1.033 | 0.903 | 0.135 | -0.006 | 0.217 |
| M003 | 277 | True | 1.058 | 1.205 | -0.131 | -0.014 | 0.728 |
| M007 | 131 | True | 1.719 | 1.08 | 0.464 | -0.029 | 0.041 |
| M013 | 60 | True | 1.851 | 1.627 | 0.129 | -0.032 | 0.378 |
| M016 | 39 | True | 1.135 | 1.358 | -0.18 | -0.106 | 0.625 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._