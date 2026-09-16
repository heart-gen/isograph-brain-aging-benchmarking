# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4476 | 0.1443 | 0.1608 | 0.84 | 0.75 | 0.94 | 1.69e-03 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2517 | 0.1331 | 0.1604 | 0.8 | 0.7 | 0.91 | 7.73e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 354 | 0.1949 | 0.153 | 1.04 | 0.78 | 1.37 | 7.97e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 2163 | 0.123 | 0.1617 | 0.76 | 0.66 | 0.88 | 1.90e-04 | logit_matched_standard |
| sQTL | all_modules | 4476 | 0.1443 | 0.1608 | 0.82 | 0.73 | 0.92 | 7.48e-04 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 2517 | 0.1331 | 0.1604 | 0.79 | 0.69 | 0.91 | 7.99e-04 | logit_matched_constraint |
| sQTL | go_invisible_modules | 354 | 0.1949 | 0.153 | 1.04 | 0.78 | 1.38 | 8.04e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 2163 | 0.123 | 0.1617 | 0.75 | 0.65 | 0.87 | 1.88e-04 | logit_matched_constraint |
| eQTL | all_modules | 5014 | 0.348 | 0.3872 | 0.84 | 0.78 | 0.9 | 1.84e-06 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2890 | 0.3343 | 0.3832 | 0.81 | 0.75 | 0.89 | 3.67e-06 | logit_matched_standard |
| eQTL | go_invisible_modules | 379 | 0.3694 | 0.3729 | 0.95 | 0.77 | 1.17 | 6.31e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 2511 | 0.329 | 0.3827 | 0.8 | 0.73 | 0.88 | 2.73e-06 | logit_matched_standard |
| eQTL | all_modules | 5014 | 0.348 | 0.3872 | 0.93 | 0.86 | 1.0 | 6.14e-02 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2890 | 0.3343 | 0.3832 | 0.91 | 0.83 | 0.99 | 3.48e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 379 | 0.3694 | 0.3729 | 0.93 | 0.75 | 1.15 | 4.93e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 2511 | 0.329 | 0.3827 | 0.91 | 0.83 | 1.0 | 5.26e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
