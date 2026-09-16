# Module-level genetic anchoring — gtex-aging/hippocampus vs GTEx Brain_Hippocampus xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region hippocampus`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 1369 | False | 0.826 | 0.713 | 0.148 | 0.001 | 0.091 |
| M004 | 272 | False | 1.346 | 1.102 | 0.201 | -0.023 | 0.148 |
| M016 | 50 | False | 0.786 | 0.842 | -0.068 | -0.071 | 0.531 |
| M019 | 40 | True | 0.856 | 1.046 | -0.2 | -0.158 | 0.592 |
| M027 | 20 | True | 1.428 | 1.391 | 0.026 | -0.885 | 0.471 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._