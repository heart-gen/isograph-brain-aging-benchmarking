# Module-level genetic anchoring — gtex-aging/cerebellum vs GTEx Brain_Cerebellum xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region cerebellum`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 349 | False | 0.901 | 0.868 | 0.037 | -0.006 | 0.408 |
| M001 | 336 | True | 0.959 | 1.165 | -0.195 | 0.0 | 0.879 |
| M002 | 234 | False | 0.846 | 0.559 | 0.413 | 0.003 | 0.027 |
| M003 | 180 | True | 0.85 | 1.286 | -0.414 | -0.022 | 0.947 |
| M009 | 103 | False | 0.616 | 0.445 | 0.326 | -0.031 | 0.121 |
| M021 | 32 | False | 0.553 | 1.083 | -0.672 | -0.033 | 0.862 |
| M024 | 29 | True | 2.053 | 1.3 | 0.457 | -0.076 | 0.173 |
| M027 | 25 | True | 0.538 | inf |  |  | NA |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._