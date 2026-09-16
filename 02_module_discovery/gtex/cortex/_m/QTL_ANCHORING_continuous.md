# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3980 | 0.2317 | 0.2219 | 0.93 | 0.9 | 0.97 | 1.19e-04 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 2551 | 0.2148 | 0.2271 | 0.9 | 0.86 | 0.94 | 1.38e-06 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 511 | 0.2329 | 0.2245 | 0.94 | 0.86 | 1.02 | 1.31e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 2040 | 0.2103 | 0.2273 | 0.9 | 0.86 | 0.94 | 8.39e-06 | ols_rankint_matched_standard |
| eQTL | all_modules | 4510 | 0.5182 | 0.5449 | 0.9 | 0.87 | 0.93 | 3.07e-09 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 2904 | 0.5017 | 0.5452 | 0.9 | 0.86 | 0.93 | 6.84e-08 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 568 | 0.5123 | 0.5391 | 0.91 | 0.84 | 0.99 | 2.90e-02 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2336 | 0.4991 | 0.544 | 0.9 | 0.86 | 0.94 | 2.03e-06 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
