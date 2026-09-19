# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5334 | 0.1149 | 0.1035 | 0.93 | 0.82 | 1.05 | 2.40e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 46 | 0.0 | 0.1094 | 0.0 | 0.0 | inf | 9.99e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 46 | 0.0 | 0.1094 | 0.0 | 0.0 | inf | 9.99e-01 | logit_matched_standard |
| sQTL | all_modules | 5334 | 0.1149 | 0.1035 | 0.89 | 0.78 | 1.01 | 7.54e-02 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 46 | 0.0 | 0.1094 | 0.0 | nan | nan | 7.95e-03 | fisher_unmatched (LinAlgError) |
| sQTL | go_visible_modules | 46 | 0.0 | 0.1094 | 0.0 | nan | nan | 7.95e-03 | fisher_unmatched (LinAlgError) |
| eQTL | all_modules | 5910 | 0.2445 | 0.2663 | 0.85 | 0.79 | 0.93 | 1.20e-04 | logit_matched_standard |
| eQTL | pheno_sig_modules | 58 | 0.1724 | 0.2572 | 0.56 | 0.28 | 1.1 | 9.32e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 58 | 0.1724 | 0.2572 | 0.56 | 0.28 | 1.1 | 9.32e-02 | logit_matched_standard |
| eQTL | all_modules | 5910 | 0.2445 | 0.2663 | 0.89 | 0.82 | 0.97 | 6.41e-03 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 58 | 0.1724 | 0.2572 | 0.59 | 0.29 | 1.17 | 1.30e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 58 | 0.1724 | 0.2572 | 0.59 | 0.29 | 1.17 | 1.30e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
