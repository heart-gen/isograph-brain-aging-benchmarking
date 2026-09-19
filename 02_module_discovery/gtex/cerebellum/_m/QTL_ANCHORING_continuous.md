# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5175 | 0.3277 | 0.2531 | 0.98 | 0.94 | 1.01 | 1.78e-01 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 2869 | 0.3283 | 0.2698 | 0.97 | 0.94 | 1.01 | 1.91e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1001 | 0.3666 | 0.2757 | 1.03 | 0.97 | 1.09 | 3.92e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 1868 | 0.3078 | 0.2785 | 0.95 | 0.91 | 0.99 | 2.82e-02 | ols_rankint_matched_standard |
| eQTL | all_modules | 5785 | 0.6145 | 0.6117 | 0.97 | 0.94 | 1.0 | 4.98e-02 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 3187 | 0.6056 | 0.614 | 0.96 | 0.93 | 1.0 | 4.59e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1098 | 0.6138 | 0.6125 | 1.01 | 0.95 | 1.07 | 8.38e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2089 | 0.6012 | 0.614 | 0.94 | 0.9 | 0.99 | 1.14e-02 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
