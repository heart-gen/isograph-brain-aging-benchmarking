# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2341 | 0.217 | 0.2257 | 0.95 | 0.91 | 1.0 | 3.79e-02 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 207 | 0.2512 | 0.2238 | 1.0 | 0.87 | 1.14 | 9.76e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 136 | 0.2868 | 0.2236 | 1.1 | 0.93 | 1.29 | 2.65e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 71 | 0.1831 | 0.2245 | 0.83 | 0.66 | 1.04 | 1.11e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 3350 | 0.463 | 0.5326 | 0.85 | 0.82 | 0.88 | 8.11e-18 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 285 | 0.407 | 0.5222 | 0.75 | 0.66 | 0.84 | 8.61e-07 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 165 | 0.5212 | 0.5205 | 0.91 | 0.78 | 1.06 | 2.17e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 120 | 0.25 | 0.5222 | 0.57 | 0.48 | 0.69 | 1.05e-09 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
