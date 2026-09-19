# Module-level genetic anchoring — gtex-aging/cortex vs GTEx Brain_Cortex xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region cortex`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 1781 | False | 0.898 | 0.955 | -0.061 | 0.0 | 0.768 |
| M002 | 862 | True | 1.021 | 1.093 | -0.068 | 0.002 | 0.749 |
| M003 | 439 | False | 0.829 | 0.849 | -0.023 | -0.002 | 0.561 |
| M021 | 26 | True | 1.32 | 1.108 | 0.176 | -0.045 | 0.380 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._