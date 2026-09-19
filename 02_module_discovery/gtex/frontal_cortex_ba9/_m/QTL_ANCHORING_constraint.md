# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6910 | 0.2123 | 0.1793 | 1.1 | 0.99 | 1.21 | 7.03e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 4741 | 0.2071 | 0.1936 | 1.0 | 0.91 | 1.11 | 9.42e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 576 | 0.2101 | 0.1987 | 0.95 | 0.77 | 1.18 | 6.47e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 4165 | 0.2067 | 0.1949 | 1.01 | 0.92 | 1.12 | 7.74e-01 | logit_matched_standard |
| sQTL | all_modules | 6910 | 0.2123 | 0.1793 | 1.07 | 0.96 | 1.19 | 2.24e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 4741 | 0.2071 | 0.1936 | 0.99 | 0.89 | 1.09 | 7.81e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 576 | 0.2101 | 0.1987 | 0.97 | 0.78 | 1.21 | 7.80e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 4165 | 0.2067 | 0.1949 | 0.99 | 0.89 | 1.1 | 8.80e-01 | logit_matched_constraint |
| eQTL | all_modules | 7700 | 0.4908 | 0.5243 | 0.86 | 0.8 | 0.92 | 1.71e-05 | logit_matched_standard |
| eQTL | pheno_sig_modules | 5304 | 0.4921 | 0.5139 | 0.92 | 0.85 | 0.98 | 1.20e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 614 | 0.4805 | 0.5066 | 0.88 | 0.75 | 1.03 | 1.19e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 4690 | 0.4936 | 0.5116 | 0.93 | 0.87 | 1.0 | 5.83e-02 | logit_matched_standard |
| eQTL | all_modules | 7700 | 0.4908 | 0.5243 | 0.91 | 0.85 | 0.98 | 1.49e-02 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 5304 | 0.4921 | 0.5139 | 0.98 | 0.91 | 1.06 | 6.34e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 614 | 0.4805 | 0.5066 | 0.87 | 0.74 | 1.03 | 1.04e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 4690 | 0.4936 | 0.5116 | 1.01 | 0.94 | 1.09 | 8.02e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
