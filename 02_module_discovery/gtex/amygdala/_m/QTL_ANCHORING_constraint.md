# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2635 | 0.1078 | 0.1092 | 0.84 | 0.73 | 0.98 | 2.35e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 976 | 0.0871 | 0.111 | 0.72 | 0.57 | 0.92 | 8.96e-03 | logit_matched_standard |
| sQTL | go_invisible_modules | 456 | 0.0965 | 0.1094 | 0.88 | 0.64 | 1.22 | 4.50e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 520 | 0.0788 | 0.1103 | 0.63 | 0.45 | 0.88 | 6.50e-03 | logit_matched_standard |
| sQTL | all_modules | 2635 | 0.1078 | 0.1092 | 0.82 | 0.7 | 0.95 | 9.69e-03 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 976 | 0.0871 | 0.111 | 0.73 | 0.57 | 0.93 | 1.14e-02 | logit_matched_constraint |
| sQTL | go_invisible_modules | 456 | 0.0965 | 0.1094 | 0.89 | 0.64 | 1.24 | 5.04e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 520 | 0.0788 | 0.1103 | 0.63 | 0.45 | 0.88 | 7.28e-03 | logit_matched_constraint |
| eQTL | all_modules | 3109 | 0.2284 | 0.2635 | 0.8 | 0.72 | 0.88 | 3.25e-06 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1198 | 0.1995 | 0.2609 | 0.69 | 0.6 | 0.8 | 1.23e-06 | logit_matched_standard |
| eQTL | go_invisible_modules | 489 | 0.2127 | 0.2572 | 0.77 | 0.61 | 0.96 | 1.85e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 709 | 0.1904 | 0.2591 | 0.66 | 0.55 | 0.8 | 3.06e-05 | logit_matched_standard |
| eQTL | all_modules | 3109 | 0.2284 | 0.2635 | 0.81 | 0.74 | 0.89 | 2.55e-05 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 1198 | 0.1995 | 0.2609 | 0.73 | 0.63 | 0.85 | 4.05e-05 | logit_matched_constraint |
| eQTL | go_invisible_modules | 489 | 0.2127 | 0.2572 | 0.82 | 0.66 | 1.03 | 9.09e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 709 | 0.1904 | 0.2591 | 0.69 | 0.57 | 0.84 | 1.63e-04 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
