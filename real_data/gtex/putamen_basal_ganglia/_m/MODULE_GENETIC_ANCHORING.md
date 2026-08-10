# Module-level genetic anchoring — gtex-aging/putamen_basal_ganglia vs GTEx Brain_Putamen_basal_ganglia xQTL

Per phenotype-associated module: the splicing-specificity contrast (log sQTL matched-OR − log eQTL matched-OR) on the module's own member genes, calibrated against a size-matched permutation null (1000 draws). A positive contrast with small `perm_p` = the module is anchored to *splicing* genetics as a unit, beyond a random equal-size gene set.

Reproduce: `python -m isograph_benchmark.real_data.module_genetic_anchoring --analysis gtex-aging --region putamen_basal_ganglia`

| module_id | module_size | go_invisible | sqtl_or | eqtl_or | contrast_log | null_mean | perm_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M009 | 256 | False | 0.837 | 0.482 | 0.552 | -0.006 | 0.006 |
| M019 | 92 | True | 1.814 | 1.083 | 0.516 | -0.036 | 0.048 |
| M027 | 61 | False | 1.164 | 1.078 | 0.076 | -0.061 | 0.390 |
| M028 | 54 | True | 2.487 | 2.861 | -0.14 | -0.039 | 0.611 |
| M029 | 51 | False | 0.834 | 0.669 | 0.22 | -0.054 | 0.291 |

_Caveat: cis-sQTL anchors member-gene splicing; a positive module contrast shows the module concentrates splicing-anchored genes above chance, not that a single variant drives the whole module. The eigenswitch×genotype and module-restricted S-LDSC tests (new controlled-genotype extraction) are the deeper follow-ons._