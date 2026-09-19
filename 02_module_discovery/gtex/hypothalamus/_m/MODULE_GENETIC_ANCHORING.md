# Module-level genetic anchoring — gtex-aging/hypothalamus vs GTEx Brain_Hypothalamus xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region hypothalamus`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M001 | 982 | True | 0.897 | 0.849 | 0.055 | 0.001 | 0.308 |
| M013 | 101 | True | 0.812 | 1.424 | -0.561 | -0.006 | 0.963 |
| M014 | 99 | False | 0.61 | 0.808 | -0.281 | -0.027 | 0.790 |
| M021 | 39 | True | 1.222 | 1.246 | -0.019 | -0.092 | 0.484 |
| M024 | 28 | True | 1.954 | 0.875 | 0.803 | -0.12 | 0.063 |
| M025 | 23 | True | 1.75 | 0.71 | 0.902 | -0.26 | 0.070 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._