# Module-level genetic anchoring — gtex-aging/hypothalamus vs GTEx Brain_Hypothalamus xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region hypothalamus`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M003 | 291 | True | 0.863 | 0.848 | 0.018 | -0.004 | 0.469 |
| M009 | 114 | True | 0.915 | 0.839 | 0.087 | -0.022 | 0.369 |
| M010 | 94 | False | 0.569 | 0.898 | -0.457 | -0.008 | 0.895 |
| M015 | 46 | True | 0.844 | 0.742 | 0.129 | -0.047 | 0.362 |
| M017 | 36 | True | 1.565 | 1.067 | 0.383 | -0.034 | 0.243 |
| M018 | 36 | True | 0.801 | 1.032 | -0.254 | -0.031 | 0.671 |
| M023 | 29 | True | 0.579 | 0.804 | -0.328 | -0.144 | 0.652 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._