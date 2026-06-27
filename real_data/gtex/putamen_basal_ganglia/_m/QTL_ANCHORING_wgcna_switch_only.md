# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Putamen_basal_ganglia xQTL (gtex-aging/putamen_basal_ganglia)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region putamen_basal_ganglia --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5848 | 0.1705 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 45 | 0.0889 | 0.1711 | 0.54 | 0.19 | 1.56 | 2.54e-01 | logit_matched |
| sQTL | go_visible_modules | 45 | 0.0889 | 0.1711 | 0.54 | 0.19 | 1.56 | 2.54e-01 | logit_matched |
| eQTL | all_modules | 7086 | 0.4515 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 93 | 0.2796 | 0.4537 | 0.45 | 0.28 | 0.7 | 5.36e-04 | logit_matched |
| eQTL | go_visible_modules | 93 | 0.2796 | 0.4537 | 0.45 | 0.28 | 0.7 | 5.36e-04 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
