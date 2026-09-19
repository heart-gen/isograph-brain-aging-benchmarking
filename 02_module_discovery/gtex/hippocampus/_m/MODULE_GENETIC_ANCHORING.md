# Module-level genetic anchoring — gtex-aging/hippocampus vs GTEx Brain_Hippocampus xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region hippocampus`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M001 | 1362 | False | 0.774 | 0.709 | 0.088 | -0.002 | 0.181 |
| M018 | 57 | True | 1.354 | 1.077 | 0.229 | -0.037 | 0.281 |
| M020 | 53 | False | 1.19 | 0.721 | 0.501 | -0.075 | 0.102 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._