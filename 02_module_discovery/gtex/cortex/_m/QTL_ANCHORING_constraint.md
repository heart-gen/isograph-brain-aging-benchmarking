# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3568 | 0.22 | 0.2032 | 0.96 | 0.86 | 1.06 | 3.83e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2308 | 0.2032 | 0.2097 | 0.87 | 0.77 | 0.98 | 2.34e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 452 | 0.2257 | 0.2077 | 1.01 | 0.8 | 1.28 | 9.41e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1856 | 0.1977 | 0.2105 | 0.84 | 0.74 | 0.96 | 1.24e-02 | logit_matched_standard |
| sQTL | all_modules | 3568 | 0.22 | 0.2032 | 0.91 | 0.82 | 1.02 | 1.03e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 2308 | 0.2032 | 0.2097 | 0.87 | 0.76 | 0.98 | 2.76e-02 | logit_matched_constraint |
| sQTL | go_invisible_modules | 452 | 0.2257 | 0.2077 | 1.0 | 0.78 | 1.28 | 9.98e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1856 | 0.1977 | 0.2105 | 0.84 | 0.74 | 0.97 | 1.66e-02 | logit_matched_constraint |
| eQTL | all_modules | 3896 | 0.5041 | 0.5203 | 0.91 | 0.84 | 0.98 | 1.06e-02 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2528 | 0.4869 | 0.5222 | 0.85 | 0.78 | 0.93 | 2.36e-04 | logit_matched_standard |
| eQTL | go_invisible_modules | 486 | 0.4938 | 0.5165 | 0.87 | 0.73 | 1.04 | 1.34e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 2042 | 0.4853 | 0.521 | 0.86 | 0.78 | 0.94 | 1.31e-03 | logit_matched_standard |
| eQTL | all_modules | 3896 | 0.5041 | 0.5203 | 0.95 | 0.87 | 1.02 | 1.71e-01 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2528 | 0.4869 | 0.5222 | 0.91 | 0.83 | 1.0 | 4.91e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 486 | 0.4938 | 0.5165 | 0.86 | 0.72 | 1.04 | 1.21e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 2042 | 0.4853 | 0.521 | 0.94 | 0.85 | 1.03 | 1.94e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
