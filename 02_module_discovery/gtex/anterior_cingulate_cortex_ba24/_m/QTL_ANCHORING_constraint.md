# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6146 | 0.1529 | 0.1559 | 0.87 | 0.78 | 0.97 | 9.66e-03 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2361 | 0.1313 | 0.1604 | 0.82 | 0.71 | 0.94 | 4.02e-03 | logit_matched_standard |
| sQTL | go_visible_modules | 2361 | 0.1313 | 0.1604 | 0.82 | 0.71 | 0.94 | 4.02e-03 | logit_matched_standard |
| sQTL | all_modules | 6146 | 0.1529 | 0.1559 | 0.84 | 0.75 | 0.94 | 2.81e-03 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 2361 | 0.1313 | 0.1604 | 0.81 | 0.7 | 0.93 | 3.74e-03 | logit_matched_constraint |
| sQTL | go_visible_modules | 2361 | 0.1313 | 0.1604 | 0.81 | 0.7 | 0.93 | 3.74e-03 | logit_matched_constraint |
| eQTL | all_modules | 6855 | 0.3559 | 0.3899 | 0.85 | 0.79 | 0.91 | 3.62e-06 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2731 | 0.3336 | 0.3826 | 0.82 | 0.75 | 0.89 | 7.76e-06 | logit_matched_standard |
| eQTL | go_visible_modules | 2731 | 0.3336 | 0.3826 | 0.82 | 0.75 | 0.89 | 7.76e-06 | logit_matched_standard |
| eQTL | all_modules | 6855 | 0.3559 | 0.3899 | 0.92 | 0.85 | 0.99 | 2.60e-02 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2731 | 0.3336 | 0.3826 | 0.91 | 0.83 | 1.0 | 6.08e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 2731 | 0.3336 | 0.3826 | 0.91 | 0.83 | 1.0 | 6.08e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
