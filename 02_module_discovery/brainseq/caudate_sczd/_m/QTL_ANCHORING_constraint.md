# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4152 | 0.2076 | 0.2091 | 1.04 | 0.95 | 1.15 | 3.92e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 378 | 0.2407 | 0.2074 | 1.11 | 0.86 | 1.43 | 4.14e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 343 | 0.242 | 0.2075 | 1.1 | 0.85 | 1.43 | 4.59e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 35 | 0.2286 | 0.2085 | 1.16 | 0.51 | 2.65 | 7.20e-01 | logit_matched_standard |
| sQTL | all_modules | 4152 | 0.2076 | 0.2091 | 1.02 | 0.92 | 1.14 | 7.16e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 378 | 0.2407 | 0.2074 | 1.02 | 0.79 | 1.31 | 9.08e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 343 | 0.242 | 0.2075 | 1.01 | 0.77 | 1.32 | 9.39e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 35 | 0.2286 | 0.2085 | 1.06 | 0.46 | 2.46 | 8.89e-01 | logit_matched_constraint |
| eQTL | all_modules | 4638 | 0.4659 | 0.526 | 0.78 | 0.72 | 0.83 | 5.47e-12 | logit_matched_standard |
| eQTL | pheno_sig_modules | 404 | 0.5322 | 0.5046 | 1.12 | 0.92 | 1.37 | 2.69e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 363 | 0.5317 | 0.5047 | 1.1 | 0.9 | 1.36 | 3.51e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 41 | 0.5366 | 0.5053 | 1.24 | 0.67 | 2.3 | 4.96e-01 | logit_matched_standard |
| eQTL | all_modules | 4638 | 0.4659 | 0.526 | 0.81 | 0.75 | 0.88 | 2.54e-07 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 404 | 0.5322 | 0.5046 | 1.13 | 0.92 | 1.38 | 2.53e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 363 | 0.5317 | 0.5047 | 1.11 | 0.89 | 1.37 | 3.54e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 41 | 0.5366 | 0.5053 | 1.3 | 0.69 | 2.44 | 4.17e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
