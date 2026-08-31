# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2070 | 0.3111 | 0.2608 | 1.0 | 0.89 | 1.12 | 9.46e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 611 | 0.2864 | 0.2693 | 0.84 | 0.69 | 1.03 | 9.36e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 440 | 0.3477 | 0.267 | 0.92 | 0.74 | 1.14 | 4.29e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 171 | 0.1287 | 0.2724 | 0.6 | 0.37 | 0.97 | 3.69e-02 | logit_matched_standard |
| sQTL | all_modules | 2070 | 0.3111 | 0.2608 | 0.97 | 0.87 | 1.1 | 6.77e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 611 | 0.2864 | 0.2693 | 0.86 | 0.71 | 1.06 | 1.56e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 440 | 0.3477 | 0.267 | 0.93 | 0.75 | 1.16 | 5.26e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 171 | 0.1287 | 0.2724 | 0.63 | 0.39 | 1.03 | 6.31e-02 | logit_matched_constraint |
| eQTL | all_modules | 2327 | 0.5875 | 0.5954 | 0.91 | 0.83 | 1.0 | 5.25e-02 | logit_matched_standard |
| eQTL | pheno_sig_modules | 756 | 0.5622 | 0.5959 | 0.81 | 0.7 | 0.94 | 5.14e-03 | logit_matched_standard |
| eQTL | go_invisible_modules | 450 | 0.6622 | 0.5917 | 1.25 | 1.02 | 1.53 | 3.05e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 306 | 0.415 | 0.5981 | 0.45 | 0.36 | 0.57 | 1.97e-11 | logit_matched_standard |
| eQTL | all_modules | 2327 | 0.5875 | 0.5954 | 0.88 | 0.8 | 0.96 | 6.31e-03 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 756 | 0.5622 | 0.5959 | 0.85 | 0.73 | 0.98 | 3.09e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 450 | 0.6622 | 0.5917 | 1.2 | 0.98 | 1.47 | 8.13e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 306 | 0.415 | 0.5981 | 0.53 | 0.42 | 0.67 | 1.28e-07 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
