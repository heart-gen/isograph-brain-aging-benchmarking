# Module-level genetic anchoring — gtex-aging/cerebellum vs GTEx Brain_Cerebellum xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region cerebellum`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 1303 | False | 0.789 | 0.907 | -0.139 | -0.002 | 0.914 |
| M002 | 316 | True | 1.149 | 1.199 | -0.042 | -0.004 | 0.583 |
| M003 | 314 | True | 0.848 | 0.887 | -0.045 | -0.006 | 0.586 |
| M006 | 225 | False | 0.869 | 0.579 | 0.406 | -0.007 | 0.025 |
| M008 | 149 | True | 1.52 | 1.36 | 0.111 | -0.019 | 0.298 |
| M010 | 135 | False | 1.185 | 1.166 | 0.016 | -0.019 | 0.452 |
| M014 | 46 | True | 2.022 | 1.391 | 0.374 | -0.036 | 0.180 |
| M016 | 27 | True | 1.42 | 1.458 | -0.026 | -0.053 | 0.498 |
| M017 | 26 | True | 0.697 | 0.979 | -0.339 | -0.066 | 0.681 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._