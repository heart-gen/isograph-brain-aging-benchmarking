# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Spinal_cord_cervical_c-1 xQTL (gtex-aging/spinal_cord_cervical_c_1)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region spinal_cord_cervical_c_1`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3203 | 0.1352 | 0.1402 | 0.86 | 0.76 | 0.97 | 1.74e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 41 | 0.122 | 0.1388 | 0.87 | 0.33 | 2.29 | 7.83e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 16 | 0.125 | 0.1388 | 1.04 | 0.22 | 4.91 | 9.58e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 25 | 0.12 | 0.1388 | 0.79 | 0.23 | 2.71 | 7.06e-01 | logit_matched_standard |
| sQTL | all_modules | 3203 | 0.1352 | 0.1402 | 0.84 | 0.74 | 0.96 | 7.92e-03 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 41 | 0.122 | 0.1388 | 0.86 | 0.32 | 2.29 | 7.64e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 16 | 0.125 | 0.1388 | 1.18 | 0.25 | 5.53 | 8.30e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 25 | 0.12 | 0.1388 | 0.72 | 0.2 | 2.52 | 6.07e-01 | logit_matched_constraint |
| eQTL | all_modules | 3877 | 0.3245 | 0.3464 | 0.85 | 0.79 | 0.93 | 1.24e-04 | logit_matched_standard |
| eQTL | pheno_sig_modules | 48 | 0.4583 | 0.3399 | 1.75 | 0.99 | 3.09 | 5.50e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 18 | 0.4444 | 0.3402 | 1.6 | 0.63 | 4.06 | 3.25e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 30 | 0.4667 | 0.34 | 1.84 | 0.9 | 3.79 | 9.62e-02 | logit_matched_standard |
| eQTL | all_modules | 3877 | 0.3245 | 0.3464 | 0.87 | 0.8 | 0.95 | 9.44e-04 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 48 | 0.4583 | 0.3399 | 1.79 | 1.0 | 3.21 | 4.95e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 18 | 0.4444 | 0.3402 | 1.82 | 0.7 | 4.73 | 2.20e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 30 | 0.4667 | 0.34 | 1.78 | 0.85 | 3.7 | 1.25e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
