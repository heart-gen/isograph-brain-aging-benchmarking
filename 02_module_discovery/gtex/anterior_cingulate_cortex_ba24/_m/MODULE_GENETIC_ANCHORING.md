# Module-level genetic anchoring — gtex-aging/anterior_cingulate_cortex_ba24 vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 2563 | False | 0.759 | 0.898 | -0.169 | -0.001 | 0.984 |
| M001 | 1458 | True | 0.958 | 0.893 | 0.07 | -0.002 | 0.239 |
| M003 | 630 | False | 1.275 | 0.988 | 0.254 | 0.002 | 0.029 |
| M007 | 234 | True | 1.228 | 0.894 | 0.317 | -0.008 | 0.068 |
| M009 | 136 | True | 0.93 | 1.19 | -0.246 | -0.011 | 0.787 |
| M019 | 25 | True | 1.313 | 1.536 | -0.157 | -0.661 | 0.557 |
| M021 | 22 | True | 0.393 | 0.874 | -0.8 | -0.447 | 0.841 |
| M023 | 21 | True | 0.726 | 1.585 | -0.781 | -0.566 | 0.821 |
| M024 | 21 | True | 1.186 | 0.677 | 0.56 | -1.154 | 0.187 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._