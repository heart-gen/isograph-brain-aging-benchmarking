# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 1971 | 0.1948 | 0.211 | 0.87 | 0.76 | 0.99 | 3.25e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 174 | 0.2414 | 0.2077 | 1.19 | 0.82 | 1.71 | 3.59e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 113 | 0.2743 | 0.2075 | 1.35 | 0.87 | 2.09 | 1.76e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 61 | 0.1803 | 0.2083 | 0.89 | 0.45 | 1.75 | 7.37e-01 | logit_matched_standard |
| sQTL | all_modules | 1971 | 0.1948 | 0.211 | 0.87 | 0.76 | 0.99 | 3.49e-02 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 174 | 0.2414 | 0.2077 | 1.19 | 0.82 | 1.72 | 3.73e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 113 | 0.2743 | 0.2075 | 1.33 | 0.85 | 2.07 | 2.11e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 61 | 0.1803 | 0.2083 | 0.92 | 0.46 | 1.83 | 8.06e-01 | logit_matched_constraint |
| eQTL | all_modules | 2296 | 0.4512 | 0.5153 | 0.75 | 0.68 | 0.82 | 3.93e-10 | logit_matched_standard |
| eQTL | pheno_sig_modules | 222 | 0.3919 | 0.5064 | 0.61 | 0.46 | 0.8 | 3.40e-04 | logit_matched_standard |
| eQTL | go_invisible_modules | 127 | 0.4882 | 0.5047 | 0.89 | 0.63 | 1.27 | 5.25e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 95 | 0.2632 | 0.5062 | 0.34 | 0.22 | 0.54 | 5.07e-06 | logit_matched_standard |
| eQTL | all_modules | 2296 | 0.4512 | 0.5153 | 0.78 | 0.71 | 0.86 | 1.94e-07 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 222 | 0.3919 | 0.5064 | 0.58 | 0.44 | 0.77 | 1.41e-04 | logit_matched_constraint |
| eQTL | go_invisible_modules | 127 | 0.4882 | 0.5047 | 0.91 | 0.63 | 1.3 | 5.94e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 95 | 0.2632 | 0.5062 | 0.31 | 0.19 | 0.49 | 6.09e-07 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
