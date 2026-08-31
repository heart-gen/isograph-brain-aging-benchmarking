# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2181 | 0.2393 | 0.201 | 1.06 | 0.94 | 1.2 | 3.14e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1238 | 0.2294 | 0.2057 | 0.99 | 0.85 | 1.15 | 9.22e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 769 | 0.2523 | 0.2051 | 1.03 | 0.86 | 1.23 | 7.38e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 469 | 0.1919 | 0.209 | 0.93 | 0.72 | 1.19 | 5.43e-01 | logit_matched_standard |
| sQTL | all_modules | 2181 | 0.2393 | 0.201 | 1.04 | 0.92 | 1.17 | 5.23e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 1238 | 0.2294 | 0.2057 | 1.02 | 0.88 | 1.19 | 7.80e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 769 | 0.2523 | 0.2051 | 1.06 | 0.88 | 1.27 | 5.49e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 469 | 0.1919 | 0.209 | 0.96 | 0.74 | 1.23 | 7.23e-01 | logit_matched_constraint |
| eQTL | all_modules | 2530 | 0.4818 | 0.5211 | 0.83 | 0.76 | 0.9 | 2.26e-05 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1416 | 0.471 | 0.5188 | 0.78 | 0.7 | 0.88 | 1.73e-05 | logit_matched_standard |
| eQTL | go_invisible_modules | 803 | 0.5019 | 0.5147 | 0.88 | 0.76 | 1.01 | 7.83e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 613 | 0.4307 | 0.5178 | 0.7 | 0.59 | 0.82 | 1.74e-05 | logit_matched_standard |
| eQTL | all_modules | 2530 | 0.4818 | 0.5211 | 0.84 | 0.77 | 0.92 | 1.28e-04 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 1416 | 0.471 | 0.5188 | 0.83 | 0.74 | 0.93 | 9.68e-04 | logit_matched_constraint |
| eQTL | go_invisible_modules | 803 | 0.5019 | 0.5147 | 0.94 | 0.81 | 1.1 | 4.55e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 613 | 0.4307 | 0.5178 | 0.71 | 0.61 | 0.84 | 7.66e-05 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
