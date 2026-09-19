# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 7881 | 0.2223 | 0.2046 | 0.96 | 0.92 | 0.99 | 9.12e-03 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 5360 | 0.2155 | 0.2144 | 0.96 | 0.92 | 0.99 | 9.54e-03 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 702 | 0.235 | 0.2138 | 1.0 | 0.92 | 1.07 | 8.97e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 4658 | 0.2125 | 0.2161 | 0.95 | 0.92 | 0.99 | 9.11e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 9201 | 0.5089 | 0.5525 | 0.89 | 0.86 | 0.92 | 9.58e-15 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 6299 | 0.5069 | 0.5427 | 0.91 | 0.88 | 0.93 | 2.65e-10 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 775 | 0.5097 | 0.5311 | 0.94 | 0.88 | 1.01 | 1.11e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 5524 | 0.5065 | 0.5407 | 0.91 | 0.88 | 0.94 | 6.03e-09 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
