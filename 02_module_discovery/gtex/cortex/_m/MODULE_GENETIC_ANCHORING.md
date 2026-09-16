# Module-level genetic anchoring — gtex-aging/cortex vs GTEx Brain_Cortex xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region cortex`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 1724 | False | 0.961 | 0.984 | -0.024 | -0.0 | 0.589 |
| M001 | 612 | False | 0.806 | 0.918 | -0.131 | -0.001 | 0.842 |
| M002 | 290 | True | 0.809 | 0.913 | -0.121 | -0.008 | 0.739 |
| M003 | 249 | True | 1.092 | 1.133 | -0.037 | 0.003 | 0.572 |
| M017 | 60 | True | 1.072 | 1.112 | -0.037 | -0.029 | 0.534 |
| M021 | 29 | True | 0.844 | 1.202 | -0.354 | -0.053 | 0.709 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._