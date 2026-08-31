# Module-level genetic anchoring — gtex-aging/hypothalamus vs GTEx Brain_Hypothalamus xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region hypothalamus`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M004 | 336 | False | 0.911 | 0.942 | -0.034 | -0.004 | 0.571 |
| M005 | 240 | True | 2.127 | 1.257 | 0.526 | -0.004 | 0.004 |
| M006 | 235 | False | 0.693 | 0.731 | -0.054 | 0.002 | 0.589 |
| M007 | 210 | True | 1.031 | 0.881 | 0.157 | -0.009 | 0.244 |
| M009 | 141 | False | 0.706 | 0.746 | -0.055 | -0.015 | 0.570 |
| M017 | 60 | True | 0.749 | 0.656 | 0.133 | -0.043 | 0.343 |
| M024 | 24 | True | 0.844 | 1.089 | -0.254 | -0.542 | 0.596 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._