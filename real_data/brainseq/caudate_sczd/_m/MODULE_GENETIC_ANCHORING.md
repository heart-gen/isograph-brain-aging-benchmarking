# Module-level genetic anchoring — brainseq-sczd vs GTEx Brain_Caudate_basal_ganglia xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis brainseq-sczd`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M010 | 83 | True | 1.663 | 1.095 | 0.418 | -0.044 | 0.139 |
| M011 | 65 | False | 0.664 | 0.844 | -0.241 | -0.045 | 0.680 |
| M012 | 58 | False | 0.821 | 0.314 | 0.96 | -0.056 | 0.019 |
| M020 | 30 | True |  | 1.088 |  |  | NA |
| M022 | 29 | True | 1.282 | 1.022 | 0.226 | -0.284 | 0.357 |
| M023 | 27 | True | 1.485 | 2.214 | -0.399 | -0.359 | 0.685 |
| M025 | 24 | True | 1.012 | 1.567 | -0.437 | -1.032 | 0.654 |
| M026 | 23 | True | 1.063 | 0.964 | 0.098 | -0.791 | 0.406 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._