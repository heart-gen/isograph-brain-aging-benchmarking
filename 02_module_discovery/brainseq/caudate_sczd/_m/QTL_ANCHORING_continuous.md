# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4607 | 0.2129 | 0.2306 | 0.96 | 0.93 | 0.99 | 2.00e-02 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 417 | 0.2518 | 0.2238 | 0.99 | 0.9 | 1.09 | 8.64e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 380 | 0.2553 | 0.2237 | 1.0 | 0.9 | 1.1 | 9.30e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 37 | 0.2162 | 0.2246 | 0.95 | 0.7 | 1.31 | 7.73e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 5816 | 0.4713 | 0.5494 | 0.84 | 0.81 | 0.86 | 5.23e-30 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 458 | 0.548 | 0.5244 | 1.01 | 0.92 | 1.11 | 7.98e-01 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 415 | 0.5494 | 0.5244 | 1.01 | 0.92 | 1.11 | 8.30e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 43 | 0.5349 | 0.525 | 1.03 | 0.76 | 1.38 | 8.66e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
