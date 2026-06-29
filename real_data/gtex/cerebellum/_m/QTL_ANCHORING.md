# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2355 | 0.3231 | 0.2735 | 0.96 | 0.87 | 1.07 | 4.92e-01 | logit_matched |
| sQTL | pheno_sig_modules | 708 | 0.2895 | 0.282 | 0.92 | 0.77 | 1.09 | 3.34e-01 | logit_matched |
| sQTL | go_invisible_modules | 462 | 0.3117 | 0.2814 | 0.93 | 0.76 | 1.16 | 5.37e-01 | logit_matched |
| sQTL | go_visible_modules | 246 | 0.248 | 0.2831 | 0.88 | 0.64 | 1.21 | 4.34e-01 | logit_matched |
| eQTL | all_modules | 2782 | 0.5945 | 0.6133 | 0.87 | 0.8 | 0.95 | 1.48e-03 | logit_matched |
| eQTL | pheno_sig_modules | 870 | 0.5264 | 0.6146 | 0.66 | 0.57 | 0.76 | 2.53e-09 | logit_matched |
| eQTL | go_invisible_modules | 495 | 0.5657 | 0.6118 | 0.79 | 0.65 | 0.94 | 8.97e-03 | logit_matched |
| eQTL | go_visible_modules | 375 | 0.4747 | 0.6133 | 0.54 | 0.44 | 0.66 | 3.74e-09 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
