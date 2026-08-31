# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4078 | 0.1577 | 0.1522 | 0.99 | 0.88 | 1.1 | 7.96e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2796 | 0.1509 | 0.1552 | 0.96 | 0.84 | 1.08 | 4.83e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1426 | 0.1711 | 0.1517 | 1.07 | 0.92 | 1.25 | 3.84e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1370 | 0.1299 | 0.1575 | 0.85 | 0.71 | 1.01 | 6.32e-02 | logit_matched_standard |
| sQTL | all_modules | 4078 | 0.1577 | 0.1522 | 0.96 | 0.86 | 1.08 | 5.21e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 2796 | 0.1509 | 0.1552 | 0.95 | 0.83 | 1.08 | 4.26e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1426 | 0.1711 | 0.1517 | 1.09 | 0.93 | 1.28 | 2.94e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1370 | 0.1299 | 0.1575 | 0.82 | 0.69 | 0.98 | 3.08e-02 | logit_matched_constraint |
| eQTL | all_modules | 4710 | 0.3348 | 0.3916 | 0.76 | 0.71 | 0.82 | 1.83e-12 | logit_matched_standard |
| eQTL | pheno_sig_modules | 3244 | 0.3194 | 0.3884 | 0.73 | 0.67 | 0.79 | 2.75e-13 | logit_matched_standard |
| eQTL | go_invisible_modules | 1569 | 0.3493 | 0.3751 | 0.88 | 0.79 | 0.99 | 2.95e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 1675 | 0.2913 | 0.3834 | 0.66 | 0.59 | 0.74 | 2.86e-13 | logit_matched_standard |
| eQTL | all_modules | 4710 | 0.3348 | 0.3916 | 0.81 | 0.75 | 0.87 | 8.08e-08 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 3244 | 0.3194 | 0.3884 | 0.78 | 0.72 | 0.86 | 5.17e-08 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1569 | 0.3493 | 0.3751 | 0.96 | 0.86 | 1.08 | 4.93e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1675 | 0.2913 | 0.3834 | 0.69 | 0.62 | 0.78 | 3.46e-10 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
