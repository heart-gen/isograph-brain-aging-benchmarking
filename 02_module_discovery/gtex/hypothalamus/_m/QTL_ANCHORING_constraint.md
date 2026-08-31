# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3526 | 0.1673 | 0.166 | 0.88 | 0.79 | 0.99 | 3.31e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 435 | 0.1218 | 0.1681 | 0.71 | 0.53 | 0.96 | 2.69e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 201 | 0.1741 | 0.1663 | 1.07 | 0.73 | 1.56 | 7.45e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 234 | 0.0769 | 0.1682 | 0.44 | 0.27 | 0.72 | 1.11e-03 | logit_matched_standard |
| sQTL | all_modules | 3526 | 0.1673 | 0.166 | 0.88 | 0.78 | 0.99 | 3.29e-02 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 435 | 0.1218 | 0.1681 | 0.72 | 0.53 | 0.98 | 3.51e-02 | logit_matched_constraint |
| sQTL | go_invisible_modules | 201 | 0.1741 | 0.1663 | 1.04 | 0.7 | 1.53 | 8.60e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 234 | 0.0769 | 0.1682 | 0.46 | 0.28 | 0.76 | 2.41e-03 | logit_matched_constraint |
| eQTL | all_modules | 4034 | 0.326 | 0.3667 | 0.8 | 0.74 | 0.87 | 3.39e-08 | logit_matched_standard |
| eQTL | pheno_sig_modules | 522 | 0.2969 | 0.3573 | 0.72 | 0.59 | 0.87 | 7.28e-04 | logit_matched_standard |
| eQTL | go_invisible_modules | 214 | 0.3224 | 0.3556 | 0.82 | 0.61 | 1.09 | 1.70e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 308 | 0.2792 | 0.3568 | 0.66 | 0.51 | 0.85 | 1.35e-03 | logit_matched_standard |
| eQTL | all_modules | 4034 | 0.326 | 0.3667 | 0.81 | 0.75 | 0.88 | 3.84e-07 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 522 | 0.2969 | 0.3573 | 0.71 | 0.59 | 0.86 | 6.06e-04 | logit_matched_constraint |
| eQTL | go_invisible_modules | 214 | 0.3224 | 0.3556 | 0.84 | 0.63 | 1.13 | 2.60e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 308 | 0.2792 | 0.3568 | 0.64 | 0.49 | 0.82 | 5.19e-04 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
