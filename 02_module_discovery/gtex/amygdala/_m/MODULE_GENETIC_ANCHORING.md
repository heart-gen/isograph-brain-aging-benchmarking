# Module-level genetic anchoring — gtex-aging/amygdala vs GTEx Brain_Amygdala xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region amygdala`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 1241 | False | 0.723 | 0.817 | -0.123 | -0.002 | 0.836 |
| M001 | 1002 | False | 0.786 | 0.977 | -0.218 | -0.004 | 0.937 |
| M003 | 449 | True | 1.004 | 1.058 | -0.052 | -0.003 | 0.601 |
| M007 | 137 | False | 2.384 | 0.85 | 1.031 | -0.035 | 0.001 |
| M008 | 115 | True | 1.147 | 0.741 | 0.437 | -0.011 | 0.096 |
| M009 | 76 | True | 1.3 | 1.145 | 0.127 | -0.06 | 0.349 |
| M010 | 72 | False | 0.446 | 0.798 | -0.581 | -0.06 | 0.877 |
| M015 | 36 | True | 1.662 | 1.095 | 0.417 | -1.074 | 0.231 |
| M016 | 31 | True | 1.966 | 0.907 | 0.774 | -0.79 | 0.091 |
| M018 | 26 | True | 0.373 | 0.17 | 0.787 | -0.928 | 0.115 |
| M022 | 20 | True | 1.856 | 1.746 | 0.061 | -3.19 | 0.463 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._