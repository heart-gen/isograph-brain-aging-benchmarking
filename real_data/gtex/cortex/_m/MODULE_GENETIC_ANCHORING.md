# Module-level genetic anchoring — gtex-aging/cortex vs GTEx Brain_Cortex xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region cortex`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M001 | 218 | False | 0.565 | 0.561 | 0.008 | -0.006 | 0.478 |
| M002 | 200 | False | 1.116 | 0.999 | 0.111 | -0.009 | 0.316 |
| M003 | 177 | False | 0.679 | 1.342 | -0.681 | -0.013 | 0.995 |
| M004 | 159 | True | 0.739 | 0.782 | -0.056 | -0.013 | 0.589 |
| M005 | 157 | True | 1.184 | 0.811 | 0.378 | -0.012 | 0.061 |
| M006 | 127 | False | 1.201 | 0.946 | 0.239 | -0.037 | 0.176 |
| M009 | 108 | True | 0.795 | 0.822 | -0.034 | -0.016 | 0.502 |
| M010 | 106 | True | 0.798 | 1.491 | -0.625 | -0.005 | 0.977 |
| M011 | 105 | True | 1.034 | 0.8 | 0.257 | -0.025 | 0.192 |
| M013 | 92 | True | 0.993 | 1.601 | -0.477 | -0.009 | 0.920 |
| M016 | 82 | True | 0.953 | 1.18 | -0.213 | -0.019 | 0.710 |
| M017 | 75 | False | 0.805 | 0.924 | -0.137 | -0.046 | 0.613 |
| M021 | 49 | True | 1.228 | 1.888 | -0.43 | -0.019 | 0.821 |
| M022 | 43 | True | 0.699 | 0.849 | -0.194 | -0.061 | 0.611 |
| M027 | 27 | True | 0.677 | 0.845 | -0.222 | -0.093 | 0.602 |
| M028 | 26 | True | 1.006 | 0.591 | 0.533 | -0.167 | 0.182 |
| M029 | 26 | True | 0.901 | 0.897 | 0.004 | -0.085 | 0.469 |
| M031 | 25 | True | 1.202 | 1.205 | -0.003 | -0.085 | 0.492 |
| M032 | 25 | True | 0.248 | 0.778 | -1.143 | -0.159 | 0.928 |
| M034 | 24 | True | 1.314 | 0.731 | 0.587 | -0.131 | 0.144 |
| M036 | 22 | True | 0.945 | 1.597 | -0.525 | -0.194 | 0.753 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._