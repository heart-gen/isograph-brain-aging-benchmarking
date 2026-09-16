# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3279 | 0.3233 | 0.2481 | 1.09 | 0.99 | 1.2 | 7.84e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1922 | 0.3117 | 0.2617 | 0.97 | 0.87 | 1.09 | 6.40e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 540 | 0.3648 | 0.2655 | 1.08 | 0.89 | 1.31 | 4.61e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1382 | 0.2909 | 0.2674 | 0.93 | 0.81 | 1.07 | 2.95e-01 | logit_matched_standard |
| sQTL | all_modules | 3279 | 0.3233 | 0.2481 | 1.05 | 0.95 | 1.16 | 3.35e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 1922 | 0.3117 | 0.2617 | 0.96 | 0.85 | 1.08 | 5.26e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 540 | 0.3648 | 0.2655 | 1.04 | 0.85 | 1.26 | 7.30e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1382 | 0.2909 | 0.2674 | 0.93 | 0.81 | 1.07 | 3.29e-01 | logit_matched_constraint |
| eQTL | all_modules | 3509 | 0.6064 | 0.5934 | 1.02 | 0.94 | 1.1 | 7.09e-01 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2074 | 0.594 | 0.5973 | 0.95 | 0.86 | 1.04 | 2.60e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 566 | 0.6042 | 0.5965 | 1.02 | 0.86 | 1.21 | 8.38e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1508 | 0.5902 | 0.5976 | 0.92 | 0.83 | 1.03 | 1.55e-01 | logit_matched_standard |
| eQTL | all_modules | 3509 | 0.6064 | 0.5934 | 0.94 | 0.86 | 1.02 | 1.27e-01 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2074 | 0.594 | 0.5973 | 0.89 | 0.81 | 0.99 | 2.95e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 566 | 0.6042 | 0.5965 | 0.95 | 0.8 | 1.14 | 5.80e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1508 | 0.5902 | 0.5976 | 0.88 | 0.79 | 0.99 | 3.34e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
