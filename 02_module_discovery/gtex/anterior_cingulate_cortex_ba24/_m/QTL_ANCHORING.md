# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4757 | 0.1682 | 0.172 | 0.87 | 0.79 | 0.96 | 5.24e-03 | logit_matched |
| sQTL | pheno_sig_modules | 3236 | 0.1582 | 0.1746 | 0.84 | 0.75 | 0.94 | 2.04e-03 | logit_matched |
| sQTL | go_invisible_modules | 1652 | 0.181 | 0.1692 | 0.96 | 0.84 | 1.11 | 5.92e-01 | logit_matched |
| sQTL | go_visible_modules | 1584 | 0.1345 | 0.1755 | 0.75 | 0.64 | 0.88 | 3.63e-04 | logit_matched |
| eQTL | all_modules | 5799 | 0.3656 | 0.426 | 0.76 | 0.71 | 0.81 | 3.45e-16 | logit_matched |
| eQTL | pheno_sig_modules | 3954 | 0.3488 | 0.4228 | 0.73 | 0.68 | 0.79 | 1.47e-16 | logit_matched |
| eQTL | go_invisible_modules | 1910 | 0.378 | 0.41 | 0.87 | 0.79 | 0.96 | 5.87e-03 | logit_matched |
| eQTL | go_visible_modules | 2044 | 0.3214 | 0.4175 | 0.67 | 0.61 | 0.74 | 9.98e-16 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
