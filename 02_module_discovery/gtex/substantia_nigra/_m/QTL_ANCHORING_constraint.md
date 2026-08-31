# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3356 | 0.1138 | 0.1086 | 0.91 | 0.8 | 1.05 | 1.96e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 48 | 0.0625 | 0.1104 | 0.51 | 0.15 | 1.69 | 2.70e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 26 | 0.1154 | 0.1102 | 0.88 | 0.25 | 3.03 | 8.35e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 22 | 0.0 | 0.1104 | 0.0 | 0.0 | inf | 1.00e+00 | logit_matched_standard |
| sQTL | all_modules | 3356 | 0.1138 | 0.1086 | 0.9 | 0.78 | 1.03 | 1.32e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 48 | 0.0625 | 0.1104 | 0.56 | 0.17 | 1.88 | 3.51e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 26 | 0.1154 | 0.1102 | 0.89 | 0.25 | 3.1 | 8.52e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 22 | 0.0 | 0.1104 | 0.0 | 0.0 | inf | 9.94e-01 | logit_matched_constraint |
| eQTL | all_modules | 3889 | 0.2376 | 0.2533 | 0.88 | 0.8 | 0.96 | 3.49e-03 | logit_matched_standard |
| eQTL | pheno_sig_modules | 75 | 0.1067 | 0.2496 | 0.35 | 0.17 | 0.72 | 4.63e-03 | logit_matched_standard |
| eQTL | go_invisible_modules | 35 | 0.0857 | 0.2493 | 0.27 | 0.08 | 0.88 | 2.93e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 40 | 0.125 | 0.2492 | 0.42 | 0.16 | 1.08 | 7.14e-02 | logit_matched_standard |
| eQTL | all_modules | 3889 | 0.2376 | 0.2533 | 0.9 | 0.82 | 0.98 | 1.74e-02 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 75 | 0.1067 | 0.2496 | 0.38 | 0.18 | 0.79 | 9.80e-03 | logit_matched_constraint |
| eQTL | go_invisible_modules | 35 | 0.0857 | 0.2493 | 0.29 | 0.09 | 0.97 | 4.40e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 40 | 0.125 | 0.2492 | 0.46 | 0.18 | 1.18 | 1.06e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
