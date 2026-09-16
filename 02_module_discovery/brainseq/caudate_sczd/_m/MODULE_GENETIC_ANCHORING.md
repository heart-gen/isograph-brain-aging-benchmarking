# Module-level genetic anchoring — brainseq-sczd vs GTEx Brain_Caudate_basal_ganglia xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis brainseq-sczd`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M010 | 145 | True | 1.151 | 1.476 | -0.249 | -0.012 | 0.789 |
| M012 | 113 | True | 0.97 | 0.739 | 0.272 | -0.028 | 0.175 |
| M013 | 105 | True | 1.256 | 1.531 | -0.198 | -0.003 | 0.715 |
| M014 | 69 | True | 0.988 | 0.764 | 0.257 | -0.036 | 0.250 |
| M015 | 66 | True | 1.034 | 1.138 | -0.096 | -0.049 | 0.560 |
| M017 | 59 | True | 1.287 | 0.693 | 0.62 | -0.045 | 0.062 |
| M018 | 53 | True | 1.779 | 1.744 | 0.02 | -0.048 | 0.458 |
| M020 | 43 | False | 0.828 | 1.234 | -0.4 | -0.024 | 0.762 |
| M030 | 22 | True | 1.935 | 0.651 | 1.089 | -0.563 | 0.046 |
| M031 | 21 | True | 0.905 | 1.528 | -0.524 | -0.366 | 0.731 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._