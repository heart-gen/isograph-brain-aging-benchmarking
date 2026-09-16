# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5196 | 0.2009 | 0.1978 | 0.99 | 0.9 | 1.09 | 8.46e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 3817 | 0.1907 | 0.2035 | 0.91 | 0.82 | 1.01 | 9.21e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 1599 | 0.2033 | 0.1986 | 0.93 | 0.81 | 1.07 | 3.31e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 2218 | 0.1817 | 0.2034 | 0.93 | 0.82 | 1.05 | 2.53e-01 | logit_matched_standard |
| sQTL | all_modules | 5196 | 0.2009 | 0.1978 | 0.96 | 0.87 | 1.07 | 4.93e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 3817 | 0.1907 | 0.2035 | 0.91 | 0.82 | 1.02 | 1.09e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1599 | 0.2033 | 0.1986 | 0.95 | 0.82 | 1.09 | 4.74e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 2218 | 0.1817 | 0.2034 | 0.92 | 0.81 | 1.05 | 2.12e-01 | logit_matched_constraint |
| eQTL | all_modules | 5805 | 0.4827 | 0.5222 | 0.85 | 0.79 | 0.91 | 3.29e-06 | logit_matched_standard |
| eQTL | pheno_sig_modules | 4302 | 0.4823 | 0.516 | 0.88 | 0.82 | 0.95 | 6.52e-04 | logit_matched_standard |
| eQTL | go_invisible_modules | 1709 | 0.4816 | 0.5088 | 0.89 | 0.8 | 0.98 | 1.97e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 2593 | 0.4828 | 0.5107 | 0.91 | 0.84 | 1.0 | 3.84e-02 | logit_matched_standard |
| eQTL | all_modules | 5805 | 0.4827 | 0.5222 | 0.91 | 0.85 | 0.98 | 1.67e-02 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 4302 | 0.4823 | 0.516 | 0.98 | 0.9 | 1.06 | 5.58e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1709 | 0.4816 | 0.5088 | 0.94 | 0.85 | 1.05 | 2.71e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 2593 | 0.4828 | 0.5107 | 1.01 | 0.92 | 1.11 | 7.90e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
