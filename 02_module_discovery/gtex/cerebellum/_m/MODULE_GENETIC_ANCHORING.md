# Module-level genetic anchoring — gtex-aging/cerebellum vs GTEx Brain_Cerebellum xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region cerebellum`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 1611 | False | 0.872 | 0.968 | -0.104 | 0.002 | 0.882 |
| M003 | 549 | True | 0.945 | 0.882 | 0.069 | -0.004 | 0.299 |
| M004 | 493 | True | 1.306 | 1.322 | -0.012 | -0.005 | 0.529 |
| M006 | 245 | False | 0.999 | 0.949 | 0.052 | -0.001 | 0.411 |
| M007 | 234 | False | 0.937 | 0.617 | 0.418 | -0.01 | 0.018 |
| M016 | 30 | True | 1.336 | 1.33 | 0.004 | -0.067 | 0.464 |
| M018 | 26 | True | 0.755 | 0.956 | -0.236 | -0.054 | 0.639 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._