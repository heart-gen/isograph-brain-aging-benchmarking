# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3741 | 0.1061 | 0.1104 | 0.84 | 0.73 | 0.96 | 8.93e-03 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2493 | 0.0951 | 0.113 | 0.77 | 0.66 | 0.9 | 1.08e-03 | logit_matched_standard |
| sQTL | go_invisible_modules | 544 | 0.1324 | 0.1078 | 0.98 | 0.75 | 1.27 | 8.60e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1949 | 0.0847 | 0.1142 | 0.73 | 0.61 | 0.87 | 3.92e-04 | logit_matched_standard |
| sQTL | all_modules | 3741 | 0.1061 | 0.1104 | 0.82 | 0.71 | 0.94 | 3.99e-03 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 2493 | 0.0951 | 0.113 | 0.75 | 0.64 | 0.88 | 4.36e-04 | logit_matched_constraint |
| sQTL | go_invisible_modules | 544 | 0.1324 | 0.1078 | 0.91 | 0.69 | 1.19 | 4.84e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1949 | 0.0847 | 0.1142 | 0.72 | 0.6 | 0.86 | 4.20e-04 | logit_matched_constraint |
| eQTL | all_modules | 4147 | 0.2274 | 0.2698 | 0.77 | 0.7 | 0.84 | 3.54e-09 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2778 | 0.2106 | 0.2687 | 0.71 | 0.64 | 0.79 | 6.62e-11 | logit_matched_standard |
| eQTL | go_invisible_modules | 583 | 0.2367 | 0.2577 | 0.86 | 0.7 | 1.04 | 1.22e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 2195 | 0.2036 | 0.2671 | 0.7 | 0.62 | 0.78 | 2.80e-10 | logit_matched_standard |
| eQTL | all_modules | 4147 | 0.2274 | 0.2698 | 0.81 | 0.74 | 0.89 | 7.69e-06 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2778 | 0.2106 | 0.2687 | 0.76 | 0.69 | 0.85 | 4.28e-07 | logit_matched_constraint |
| eQTL | go_invisible_modules | 583 | 0.2367 | 0.2577 | 0.82 | 0.67 | 0.99 | 4.39e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 2195 | 0.2036 | 0.2671 | 0.77 | 0.68 | 0.86 | 9.42e-06 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
