# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5021 | 0.2173 | 0.2016 | 1.11 | 1.01 | 1.22 | 3.47e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1393 | 0.2685 | 0.2002 | 1.34 | 1.17 | 1.53 | 1.73e-05 | logit_matched_standard |
| sQTL | go_invisible_modules | 1108 | 0.2581 | 0.2032 | 1.19 | 1.02 | 1.38 | 2.31e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 285 | 0.3088 | 0.206 | 1.94 | 1.48 | 2.55 | 1.51e-06 | logit_matched_standard |
| sQTL | all_modules | 5021 | 0.2173 | 0.2016 | 1.09 | 0.98 | 1.21 | 1.05e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 1393 | 0.2685 | 0.2002 | 1.31 | 1.14 | 1.5 | 1.12e-04 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1108 | 0.2581 | 0.2032 | 1.15 | 0.99 | 1.34 | 7.36e-02 | logit_matched_constraint |
| sQTL | go_visible_modules | 285 | 0.3088 | 0.206 | 1.98 | 1.5 | 2.62 | 1.39e-06 | logit_matched_constraint |
| eQTL | all_modules | 5579 | 0.4793 | 0.5237 | 0.82 | 0.77 | 0.88 | 4.27e-08 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1509 | 0.5156 | 0.5041 | 1.04 | 0.94 | 1.16 | 4.34e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 1176 | 0.523 | 0.5037 | 1.07 | 0.94 | 1.2 | 2.99e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 333 | 0.4895 | 0.5058 | 0.97 | 0.78 | 1.2 | 7.64e-01 | logit_matched_standard |
| eQTL | all_modules | 5579 | 0.4793 | 0.5237 | 0.87 | 0.8 | 0.93 | 2.20e-04 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 1509 | 0.5156 | 0.5041 | 1.05 | 0.94 | 1.18 | 3.55e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1176 | 0.523 | 0.5037 | 1.06 | 0.94 | 1.2 | 3.38e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 333 | 0.4895 | 0.5058 | 1.01 | 0.81 | 1.27 | 9.00e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
