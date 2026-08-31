# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4757 | 0.1682 | 0.172 | 0.97 | 0.93 | 1.0 | 4.62e-02 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 3236 | 0.1582 | 0.1746 | 0.95 | 0.91 | 0.99 | 8.84e-03 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1652 | 0.181 | 0.1692 | 0.97 | 0.92 | 1.02 | 1.93e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 1584 | 0.1345 | 0.1755 | 0.95 | 0.9 | 1.0 | 3.34e-02 | ols_rankint_matched_standard |
| eQTL | all_modules | 5799 | 0.3656 | 0.426 | 0.85 | 0.82 | 0.88 | 8.00e-24 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 3954 | 0.3488 | 0.4228 | 0.84 | 0.81 | 0.87 | 6.68e-22 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1910 | 0.378 | 0.41 | 0.92 | 0.88 | 0.96 | 6.01e-04 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2044 | 0.3214 | 0.4175 | 0.81 | 0.77 | 0.84 | 5.94e-20 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
