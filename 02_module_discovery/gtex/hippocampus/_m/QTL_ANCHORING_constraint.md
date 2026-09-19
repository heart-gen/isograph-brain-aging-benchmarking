# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5871 | 0.1524 | 0.1477 | 0.84 | 0.75 | 0.94 | 1.79e-03 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1053 | 0.1225 | 0.153 | 0.76 | 0.62 | 0.93 | 6.73e-03 | logit_matched_standard |
| sQTL | go_invisible_modules | 32 | 0.1875 | 0.1501 | 1.12 | 0.45 | 2.8 | 8.14e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1021 | 0.1205 | 0.1531 | 0.75 | 0.61 | 0.92 | 5.13e-03 | logit_matched_standard |
| sQTL | all_modules | 5871 | 0.1524 | 0.1477 | 0.82 | 0.73 | 0.92 | 5.57e-04 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 1053 | 0.1225 | 0.153 | 0.76 | 0.62 | 0.93 | 7.58e-03 | logit_matched_constraint |
| sQTL | go_invisible_modules | 32 | 0.1875 | 0.1501 | 0.96 | 0.38 | 2.43 | 9.26e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1021 | 0.1205 | 0.1531 | 0.75 | 0.61 | 0.92 | 7.00e-03 | logit_matched_constraint |
| eQTL | all_modules | 6497 | 0.3522 | 0.3723 | 0.87 | 0.81 | 0.93 | 1.43e-04 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1237 | 0.2716 | 0.3718 | 0.64 | 0.56 | 0.73 | 1.61e-11 | logit_matched_standard |
| eQTL | go_invisible_modules | 32 | 0.375 | 0.3627 | 1.05 | 0.51 | 2.15 | 8.98e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1205 | 0.2689 | 0.3718 | 0.63 | 0.55 | 0.72 | 8.07e-12 | logit_matched_standard |
| eQTL | all_modules | 6497 | 0.3522 | 0.3723 | 0.89 | 0.83 | 0.96 | 3.13e-03 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 1237 | 0.2716 | 0.3718 | 0.69 | 0.6 | 0.79 | 8.31e-08 | logit_matched_constraint |
| eQTL | go_invisible_modules | 32 | 0.375 | 0.3627 | 0.98 | 0.48 | 2.04 | 9.66e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1205 | 0.2689 | 0.3718 | 0.69 | 0.6 | 0.79 | 6.04e-08 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
