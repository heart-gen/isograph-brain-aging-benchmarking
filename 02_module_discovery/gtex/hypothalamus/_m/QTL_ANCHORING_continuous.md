# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6581 | 0.1854 | 0.1796 | 0.94 | 0.91 | 0.97 | 7.55e-05 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 1111 | 0.18 | 0.1826 | 0.93 | 0.88 | 0.99 | 2.72e-02 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1039 | 0.1829 | 0.1823 | 0.94 | 0.88 | 1.0 | 4.02e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 72 | 0.1389 | 0.1826 | 0.91 | 0.73 | 1.14 | 4.09e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 7457 | 0.3731 | 0.3976 | 0.93 | 0.91 | 0.96 | 7.95e-06 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 1258 | 0.3498 | 0.3904 | 0.92 | 0.86 | 0.97 | 2.51e-03 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1159 | 0.3529 | 0.39 | 0.91 | 0.86 | 0.97 | 3.53e-03 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 99 | 0.3131 | 0.3881 | 0.93 | 0.76 | 1.13 | 4.58e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
