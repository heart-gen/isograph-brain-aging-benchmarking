# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5145 | 0.2278 | 0.1926 | 1.01 | 0.92 | 1.11 | 8.47e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2417 | 0.2168 | 0.2062 | 0.94 | 0.83 | 1.05 | 2.68e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 717 | 0.2594 | 0.205 | 1.08 | 0.9 | 1.3 | 4.18e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1700 | 0.1988 | 0.2101 | 0.88 | 0.77 | 1.01 | 6.36e-02 | logit_matched_standard |
| sQTL | all_modules | 5145 | 0.2278 | 0.1926 | 0.96 | 0.87 | 1.06 | 4.21e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 2417 | 0.2168 | 0.2062 | 0.92 | 0.82 | 1.04 | 2.07e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 717 | 0.2594 | 0.205 | 1.06 | 0.88 | 1.28 | 5.36e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1700 | 0.1988 | 0.2101 | 0.87 | 0.76 | 1.0 | 5.70e-02 | logit_matched_constraint |
| eQTL | all_modules | 5610 | 0.5119 | 0.5183 | 0.94 | 0.87 | 1.0 | 6.47e-02 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2658 | 0.4929 | 0.5212 | 0.87 | 0.8 | 0.95 | 1.17e-03 | logit_matched_standard |
| eQTL | go_invisible_modules | 749 | 0.5407 | 0.5143 | 1.05 | 0.9 | 1.22 | 5.33e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1909 | 0.4741 | 0.5224 | 0.82 | 0.74 | 0.9 | 4.26e-05 | logit_matched_standard |
| eQTL | all_modules | 5610 | 0.5119 | 0.5183 | 0.97 | 0.9 | 1.04 | 3.81e-01 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 2658 | 0.4929 | 0.5212 | 0.91 | 0.84 | 1.0 | 4.92e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 749 | 0.5407 | 0.5143 | 1.06 | 0.91 | 1.23 | 4.87e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1909 | 0.4741 | 0.5224 | 0.87 | 0.79 | 0.96 | 7.02e-03 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
