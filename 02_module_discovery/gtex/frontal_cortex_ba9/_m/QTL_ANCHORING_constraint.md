# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5931 | 0.2054 | 0.1923 | 1.02 | 0.93 | 1.13 | 6.20e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 5110 | 0.2029 | 0.1959 | 1.0 | 0.91 | 1.1 | 9.67e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 2206 | 0.2335 | 0.1909 | 1.11 | 0.99 | 1.25 | 7.21e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 2904 | 0.1798 | 0.2056 | 0.91 | 0.81 | 1.01 | 8.45e-02 | logit_matched_standard |
| sQTL | all_modules | 5931 | 0.2054 | 0.1923 | 1.03 | 0.93 | 1.14 | 6.19e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 5110 | 0.2029 | 0.1959 | 1.02 | 0.92 | 1.13 | 7.35e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 2206 | 0.2335 | 0.1909 | 1.14 | 1.01 | 1.29 | 3.09e-02 | logit_matched_constraint |
| sQTL | go_visible_modules | 2904 | 0.1798 | 0.2056 | 0.91 | 0.81 | 1.02 | 9.69e-02 | logit_matched_constraint |
| eQTL | all_modules | 6761 | 0.4792 | 0.5283 | 0.81 | 0.75 | 0.86 | 5.42e-10 | logit_matched_standard |
| eQTL | pheno_sig_modules | 5842 | 0.4812 | 0.5211 | 0.84 | 0.79 | 0.9 | 9.16e-07 | logit_matched_standard |
| eQTL | go_invisible_modules | 2445 | 0.5072 | 0.5035 | 1.0 | 0.92 | 1.09 | 9.76e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 3397 | 0.4625 | 0.5178 | 0.8 | 0.74 | 0.86 | 1.82e-08 | logit_matched_standard |
| eQTL | all_modules | 6761 | 0.4792 | 0.5283 | 0.86 | 0.8 | 0.92 | 3.73e-05 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 5842 | 0.4812 | 0.5211 | 0.91 | 0.84 | 0.97 | 6.71e-03 | logit_matched_constraint |
| eQTL | go_invisible_modules | 2445 | 0.5072 | 0.5035 | 1.03 | 0.94 | 1.13 | 4.71e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 3397 | 0.4625 | 0.5178 | 0.86 | 0.79 | 0.93 | 2.07e-04 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
