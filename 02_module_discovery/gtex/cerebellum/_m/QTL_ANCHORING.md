# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2355 | 0.3231 | 0.2735 | 0.96 | 0.87 | 1.07 | 4.92e-01 | logit_matched |
| sQTL | pheno_sig_modules | 700 | 0.3071 | 0.281 | 0.86 | 0.72 | 1.03 | 1.07e-01 | logit_matched |
| sQTL | go_invisible_modules | 504 | 0.369 | 0.279 | 0.94 | 0.77 | 1.15 | 5.66e-01 | logit_matched |
| sQTL | go_visible_modules | 196 | 0.148 | 0.2845 | 0.61 | 0.4 | 0.93 | 2.06e-02 | logit_matched |
| eQTL | all_modules | 2782 | 0.5945 | 0.6133 | 0.87 | 0.8 | 0.95 | 1.48e-03 | logit_matched |
| eQTL | pheno_sig_modules | 921 | 0.5765 | 0.6123 | 0.81 | 0.7 | 0.92 | 1.75e-03 | logit_matched |
| eQTL | go_invisible_modules | 552 | 0.6739 | 0.6086 | 1.24 | 1.04 | 1.49 | 1.93e-02 | logit_matched |
| eQTL | go_visible_modules | 369 | 0.4309 | 0.6141 | 0.45 | 0.37 | 0.56 | 7.79e-14 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
