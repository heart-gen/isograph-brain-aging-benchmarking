# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4233 | 0.1857 | 0.1808 | 0.93 | 0.9 | 0.96 | 6.33e-05 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 482 | 0.1349 | 0.184 | 0.85 | 0.78 | 0.93 | 2.86e-04 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 216 | 0.1806 | 0.1823 | 0.86 | 0.75 | 0.98 | 2.21e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 266 | 0.0977 | 0.1839 | 0.85 | 0.75 | 0.95 | 5.55e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 5183 | 0.3556 | 0.3985 | 0.9 | 0.87 | 0.93 | 6.03e-10 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 608 | 0.3026 | 0.3895 | 0.82 | 0.76 | 0.89 | 2.35e-06 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 234 | 0.3333 | 0.3874 | 0.88 | 0.78 | 1.0 | 5.83e-02 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 374 | 0.2834 | 0.3888 | 0.79 | 0.72 | 0.88 | 7.85e-06 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
