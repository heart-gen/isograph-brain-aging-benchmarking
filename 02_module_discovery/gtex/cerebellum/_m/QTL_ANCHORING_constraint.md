# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4614 | 0.3208 | 0.2343 | 1.11 | 1.01 | 1.21 | 3.38e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2613 | 0.3222 | 0.2543 | 1.06 | 0.96 | 1.18 | 2.63e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 883 | 0.359 | 0.2627 | 1.11 | 0.95 | 1.3 | 1.75e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1730 | 0.3035 | 0.2642 | 1.02 | 0.9 | 1.15 | 7.91e-01 | logit_matched_standard |
| sQTL | all_modules | 4614 | 0.3208 | 0.2343 | 1.06 | 0.97 | 1.17 | 2.12e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 2613 | 0.3222 | 0.2543 | 1.03 | 0.93 | 1.15 | 5.44e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 883 | 0.359 | 0.2627 | 1.07 | 0.91 | 1.25 | 4.34e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1730 | 0.3035 | 0.2642 | 1.01 | 0.89 | 1.14 | 9.22e-01 | logit_matched_constraint |
| eQTL | all_modules | 4946 | 0.6108 | 0.5886 | 1.06 | 0.98 | 1.14 | 1.43e-01 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2810 | 0.6014 | 0.5956 | 0.99 | 0.91 | 1.08 | 8.90e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 930 | 0.6129 | 0.5956 | 1.07 | 0.93 | 1.23 | 3.15e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1880 | 0.5957 | 0.597 | 0.96 | 0.86 | 1.06 | 3.69e-01 | logit_matched_standard |
| eQTL | all_modules | 4946 | 0.6108 | 0.5886 | 0.97 | 0.9 | 1.05 | 4.22e-01 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2810 | 0.6014 | 0.5956 | 0.93 | 0.85 | 1.01 | 8.93e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 930 | 0.6129 | 0.5956 | 1.0 | 0.87 | 1.15 | 9.75e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1880 | 0.5957 | 0.597 | 0.9 | 0.81 | 1.0 | 5.07e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
