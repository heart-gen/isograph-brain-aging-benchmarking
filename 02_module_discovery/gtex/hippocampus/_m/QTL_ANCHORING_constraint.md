# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4502 | 0.1519 | 0.1488 | 0.89 | 0.8 | 1.0 | 4.72e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1906 | 0.1432 | 0.1514 | 0.85 | 0.73 | 0.98 | 2.64e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 798 | 0.1667 | 0.1488 | 0.96 | 0.78 | 1.17 | 6.82e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1108 | 0.1264 | 0.1526 | 0.78 | 0.65 | 0.95 | 1.36e-02 | logit_matched_standard |
| sQTL | all_modules | 4502 | 0.1519 | 0.1488 | 0.88 | 0.79 | 0.99 | 3.37e-02 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 1906 | 0.1432 | 0.1514 | 0.83 | 0.72 | 0.97 | 1.70e-02 | logit_matched_constraint |
| sQTL | go_invisible_modules | 798 | 0.1667 | 0.1488 | 0.93 | 0.76 | 1.14 | 4.81e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1108 | 0.1264 | 0.1526 | 0.79 | 0.65 | 0.96 | 1.60e-02 | logit_matched_constraint |
| eQTL | all_modules | 5259 | 0.3282 | 0.3811 | 0.75 | 0.7 | 0.81 | 2.96e-14 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2294 | 0.296 | 0.3739 | 0.68 | 0.62 | 0.75 | 2.50e-14 | logit_matched_standard |
| eQTL | go_invisible_modules | 870 | 0.3425 | 0.3623 | 0.88 | 0.76 | 1.01 | 7.66e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 1424 | 0.2676 | 0.3717 | 0.61 | 0.54 | 0.69 | 3.32e-15 | logit_matched_standard |
| eQTL | all_modules | 5259 | 0.3282 | 0.3811 | 0.76 | 0.71 | 0.82 | 1.34e-12 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2294 | 0.296 | 0.3739 | 0.69 | 0.63 | 0.77 | 7.69e-13 | logit_matched_constraint |
| eQTL | go_invisible_modules | 870 | 0.3425 | 0.3623 | 0.88 | 0.76 | 1.02 | 8.65e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 1424 | 0.2676 | 0.3717 | 0.62 | 0.55 | 0.71 | 1.82e-13 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
