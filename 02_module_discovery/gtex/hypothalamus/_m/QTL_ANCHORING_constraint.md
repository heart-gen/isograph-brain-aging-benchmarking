# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5744 | 0.1748 | 0.1587 | 0.94 | 0.85 | 1.04 | 2.39e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 983 | 0.1597 | 0.1672 | 0.88 | 0.73 | 1.05 | 1.62e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 921 | 0.1607 | 0.1671 | 0.88 | 0.73 | 1.07 | 2.06e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 62 | 0.1452 | 0.1667 | 0.78 | 0.37 | 1.63 | 5.10e-01 | logit_matched_standard |
| sQTL | all_modules | 5744 | 0.1748 | 0.1587 | 0.91 | 0.82 | 1.01 | 8.52e-02 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 983 | 0.1597 | 0.1672 | 0.9 | 0.74 | 1.08 | 2.58e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 921 | 0.1607 | 0.1671 | 0.9 | 0.74 | 1.09 | 2.78e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 62 | 0.1452 | 0.1667 | 0.89 | 0.43 | 1.88 | 7.67e-01 | logit_matched_constraint |
| eQTL | all_modules | 6226 | 0.3516 | 0.3616 | 0.92 | 0.85 | 0.99 | 2.02e-02 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1079 | 0.3346 | 0.359 | 0.86 | 0.76 | 0.99 | 3.10e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 1001 | 0.3397 | 0.3585 | 0.88 | 0.77 | 1.01 | 7.07e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 78 | 0.2692 | 0.3576 | 0.69 | 0.41 | 1.13 | 1.42e-01 | logit_matched_standard |
| eQTL | all_modules | 6226 | 0.3516 | 0.3616 | 0.95 | 0.88 | 1.02 | 1.50e-01 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 1079 | 0.3346 | 0.359 | 0.88 | 0.77 | 1.01 | 7.11e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1001 | 0.3397 | 0.3585 | 0.91 | 0.79 | 1.05 | 2.04e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 78 | 0.2692 | 0.3576 | 0.59 | 0.36 | 0.99 | 4.39e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
