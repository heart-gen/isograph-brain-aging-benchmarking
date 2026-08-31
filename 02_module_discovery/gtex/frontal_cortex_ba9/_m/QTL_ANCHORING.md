# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6830 | 0.2151 | 0.2145 | 0.91 | 0.83 | 0.99 | 2.79e-02 | logit_matched |
| sQTL | pheno_sig_modules | 5822 | 0.2133 | 0.2159 | 0.9 | 0.83 | 0.99 | 2.13e-02 | logit_matched |
| sQTL | go_invisible_modules | 2615 | 0.2528 | 0.2058 | 1.1 | 0.99 | 1.23 | 6.58e-02 | logit_matched |
| sQTL | go_visible_modules | 3207 | 0.1812 | 0.2251 | 0.79 | 0.71 | 0.88 | 8.77e-06 | logit_matched |
| eQTL | all_modules | 8211 | 0.5007 | 0.5542 | 0.79 | 0.74 | 0.83 | 1.99e-15 | logit_matched |
| eQTL | pheno_sig_modules | 6983 | 0.5019 | 0.5475 | 0.82 | 0.77 | 0.87 | 1.25e-10 | logit_matched |
| eQTL | go_invisible_modules | 3138 | 0.5376 | 0.5285 | 1.01 | 0.94 | 1.1 | 7.38e-01 | logit_matched |
| eQTL | go_visible_modules | 3845 | 0.4728 | 0.5454 | 0.75 | 0.7 | 0.8 | 2.37e-15 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
