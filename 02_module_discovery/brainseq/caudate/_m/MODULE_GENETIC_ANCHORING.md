# Module-level genetic anchoring — brainseq-aging/caudate vs GTEx Brain_Caudate_basal_ganglia xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis brainseq-aging --region caudate`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M002 | 561 | True | 1.41 | 1.092 | 0.255 | -0.004 | 0.036 |
| M003 | 495 | True | 1.009 | 0.981 | 0.028 | -0.008 | 0.424 |
| M005 | 331 | True | 1.637 | 0.956 | 0.538 | -0.01 | 0.002 |
| M010 | 159 | True | 0.509 | 0.941 | -0.615 | -0.024 | 0.980 |
| M012 | 131 | True | 0.698 | 0.938 | -0.295 | -0.01 | 0.831 |
| M014 | 114 | True | 1.407 | 1.105 | 0.241 | -0.014 | 0.221 |
| M015 | 90 | True | 0.827 | 1.125 | -0.308 | -0.044 | 0.756 |
| M018 | 70 | True | 1.154 | 0.741 | 0.443 | -0.046 | 0.124 |
| M019 | 59 | True | 1.338 | 1.766 | -0.277 | -0.028 | 0.726 |
| M026 | 49 | True | 0.761 | 1.627 | -0.76 | -0.031 | 0.920 |
| M033 | 29 | True | 1.12 | 1.604 | -0.359 | -0.77 | 0.679 |
| M035 | 27 | True | 2.083 | 1.237 | 0.521 | -0.227 | 0.197 |
| M041 | 23 | True | 0.319 | 1.724 | -1.689 | -1.375 | 0.957 |
| M042 | 23 | True | 1.756 | 0.646 | 1.0 | -0.461 | 0.064 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._