# Module-level genetic anchoring — gtex-aging/hippocampus vs GTEx Brain_Hippocampus xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region hippocampus`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M001 | 779 | False | 1.117 | 0.912 | 0.203 | -0.003 | 0.054 |
| M002 | 494 | False | 0.707 | 0.736 | -0.04 | -0.007 | 0.575 |
| M004 | 462 | True | 1.052 | 1.108 | -0.051 | -0.008 | 0.625 |
| M005 | 430 | True | 1.268 | 1.049 | 0.19 | -0.004 | 0.120 |
| M009 | 220 | False | 0.726 | 0.81 | -0.109 | -0.024 | 0.650 |
| M010 | 209 | False | 1.055 | 0.866 | 0.197 | -0.008 | 0.203 |
| M012 | 159 | False | 0.83 | 0.717 | 0.147 | -0.019 | 0.280 |
| M021 | 65 | True | 1.139 | 1.697 | -0.399 | -0.018 | 0.836 |
| M024 | 47 | True | 0.332 | 1.025 | -1.128 | -0.032 | 0.973 |
| M026 | 43 | False | 1.461 | 0.332 | 1.483 | -0.066 | 0.001 |
| M028 | 41 | True | 1.34 | 0.908 | 0.39 | -0.069 | 0.203 |
| M029 | 35 | False | 2.105 | 0.89 | 0.861 | -0.242 | 0.046 |
| M033 | 29 | True | 0.563 | 1.336 | -0.865 | -0.392 | 0.880 |
| M038 | 20 | True | 1.448 | 1.006 | 0.364 | -1.09 | 0.275 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._