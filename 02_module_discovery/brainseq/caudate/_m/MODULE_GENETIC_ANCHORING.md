# Module-level genetic anchoring — brainseq-aging/caudate vs GTEx Brain_Caudate_basal_ganglia xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis brainseq-aging --region caudate`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M000 | 797 | False | 0.952 | 0.989 | -0.039 | -0.0 | 0.628 |
| M005 | 433 | False | 1.114 | 1.138 | -0.022 | -0.001 | 0.550 |
| M008 | 205 | True | 1.105 | 0.881 | 0.226 | -0.013 | 0.172 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._