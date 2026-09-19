# Module-level genetic anchoring — brainseq-sczd vs GTEx Brain_Caudate_basal_ganglia xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis brainseq-sczd`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M004 | 432 | True | 1.111 | 0.943 | 0.164 | -0.005 | 0.148 |
| M005 | 382 | False | 1.522 | 0.963 | 0.458 | -0.016 | 0.006 |
| M006 | 343 | True | 1.143 | 1.296 | -0.125 | 0.006 | 0.763 |
| M007 | 266 | True | 1.26 | 1.467 | -0.152 | 0.003 | 0.773 |
| M008 | 230 | True | 1.139 | 1.269 | -0.108 | -0.013 | 0.671 |
| M009 | 170 | True | 1.348 | 1.17 | 0.142 | -0.019 | 0.268 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._