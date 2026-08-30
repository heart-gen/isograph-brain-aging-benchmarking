# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2527 | 0.2489 | 0.2196 | 0.98 | 0.88 | 1.09 | 6.88e-01 | logit_matched |
| sQTL | pheno_sig_modules | 1212 | 0.2467 | 0.2229 | 0.98 | 0.85 | 1.13 | 7.68e-01 | logit_matched |
| sQTL | go_invisible_modules | 963 | 0.2451 | 0.2235 | 0.95 | 0.8 | 1.11 | 5.02e-01 | logit_matched |
| sQTL | go_visible_modules | 249 | 0.253 | 0.2245 | 1.11 | 0.82 | 1.5 | 5.05e-01 | logit_matched |
| eQTL | all_modules | 3164 | 0.5035 | 0.5442 | 0.82 | 0.76 | 0.89 | 4.66e-07 | logit_matched |
| eQTL | pheno_sig_modules | 1514 | 0.4993 | 0.5406 | 0.82 | 0.74 | 0.91 | 2.33e-04 | logit_matched |
| eQTL | go_invisible_modules | 1198 | 0.4917 | 0.5404 | 0.8 | 0.71 | 0.9 | 2.37e-04 | logit_matched |
| eQTL | go_visible_modules | 316 | 0.5285 | 0.5373 | 0.91 | 0.73 | 1.14 | 4.29e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
