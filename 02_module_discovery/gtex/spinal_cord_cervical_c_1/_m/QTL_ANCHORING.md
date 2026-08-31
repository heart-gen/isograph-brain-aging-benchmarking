# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Spinal_cord_cervical_c-1 xQTL (gtex-aging/spinal_cord_cervical_c_1)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region spinal_cord_cervical_c_1`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3698 | 0.146 | 0.1574 | 0.76 | 0.68 | 0.85 | 1.41e-06 | logit_matched |
| sQTL | pheno_sig_modules | 46 | 0.1304 | 0.1543 | 0.84 | 0.35 | 2.04 | 6.98e-01 | logit_matched |
| sQTL | go_invisible_modules | 18 | 0.1667 | 0.1542 | 1.16 | 0.32 | 4.28 | 8.20e-01 | logit_matched |
| sQTL | go_visible_modules | 28 | 0.1071 | 0.1543 | 0.66 | 0.19 | 2.24 | 5.04e-01 | logit_matched |
| eQTL | all_modules | 4761 | 0.3398 | 0.3789 | 0.79 | 0.74 | 0.85 | 8.66e-11 | logit_matched |
| eQTL | pheno_sig_modules | 59 | 0.4407 | 0.3685 | 1.41 | 0.84 | 2.36 | 1.93e-01 | logit_matched |
| eQTL | go_invisible_modules | 21 | 0.5238 | 0.3686 | 1.92 | 0.81 | 4.53 | 1.37e-01 | logit_matched |
| eQTL | go_visible_modules | 38 | 0.3947 | 0.3687 | 1.18 | 0.61 | 2.27 | 6.18e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
