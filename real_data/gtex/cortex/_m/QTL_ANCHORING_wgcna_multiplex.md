# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13376 | 0.2231 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7243 | 0.2369 | 0.2068 | 1.17 | 1.07 | 1.27 | 5.13e-04 | logit_matched |
| sQTL | go_invisible_modules | 1431 | 0.2739 | 0.217 | 1.08 | 0.95 | 1.23 | 2.54e-01 | logit_matched |
| sQTL | go_visible_modules | 6524 | 0.2333 | 0.2134 | 1.18 | 1.08 | 1.28 | 1.97e-04 | logit_matched |
| eQTL | all_modules | 17847 | 0.5335 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 9134 | 0.5504 | 0.5159 | 1.15 | 1.08 | 1.22 | 5.38e-06 | logit_matched |
| eQTL | go_invisible_modules | 1724 | 0.5505 | 0.5317 | 1.04 | 0.94 | 1.15 | 4.91e-01 | logit_matched |
| eQTL | go_visible_modules | 8236 | 0.5501 | 0.5193 | 1.14 | 1.08 | 1.21 | 1.24e-05 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
