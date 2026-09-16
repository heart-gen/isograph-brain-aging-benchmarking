# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4112 | 0.1449 | 0.1532 | 0.82 | 0.73 | 0.92 | 5.47e-04 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1263 | 0.1338 | 0.1523 | 0.8 | 0.67 | 0.95 | 1.21e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 32 | 0.0938 | 0.1503 | 0.82 | 0.24 | 2.75 | 7.42e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1231 | 0.1348 | 0.1521 | 0.8 | 0.67 | 0.95 | 1.31e-02 | logit_matched_standard |
| sQTL | all_modules | 4112 | 0.1449 | 0.1532 | 0.81 | 0.72 | 0.91 | 3.96e-04 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 1263 | 0.1338 | 0.1523 | 0.79 | 0.66 | 0.95 | 1.13e-02 | logit_matched_constraint |
| sQTL | go_invisible_modules | 32 | 0.0938 | 0.1503 | 0.96 | 0.28 | 3.25 | 9.50e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1231 | 0.1348 | 0.1521 | 0.79 | 0.66 | 0.95 | 1.09e-02 | logit_matched_constraint |
| eQTL | all_modules | 4560 | 0.3386 | 0.3748 | 0.82 | 0.76 | 0.88 | 3.26e-07 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1462 | 0.2914 | 0.3713 | 0.7 | 0.62 | 0.79 | 2.89e-09 | logit_matched_standard |
| eQTL | go_invisible_modules | 37 | 0.2973 | 0.3629 | 0.71 | 0.35 | 1.43 | 3.33e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1425 | 0.2912 | 0.3711 | 0.7 | 0.62 | 0.79 | 5.21e-09 | logit_matched_standard |
| eQTL | all_modules | 4560 | 0.3386 | 0.3748 | 0.86 | 0.79 | 0.93 | 1.85e-04 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 1462 | 0.2914 | 0.3713 | 0.75 | 0.66 | 0.84 | 2.60e-06 | logit_matched_constraint |
| eQTL | go_invisible_modules | 37 | 0.2973 | 0.3629 | 0.76 | 0.37 | 1.57 | 4.63e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1425 | 0.2912 | 0.3711 | 0.75 | 0.66 | 0.84 | 3.65e-06 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
