# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3886 | 0.1225 | 0.1252 | 0.95 | 0.91 | 0.98 | 4.29e-03 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 58 | 0.0862 | 0.1246 | 1.0 | 0.78 | 1.29 | 9.93e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 30 | 0.1 | 0.1245 | 0.92 | 0.65 | 1.31 | 6.42e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 28 | 0.0714 | 0.1245 | 1.1 | 0.76 | 1.58 | 6.21e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 4751 | 0.2572 | 0.2896 | 0.92 | 0.89 | 0.95 | 2.30e-06 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 98 | 0.1633 | 0.2818 | 0.81 | 0.67 | 0.99 | 4.27e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 40 | 0.1 | 0.2815 | 0.83 | 0.61 | 1.13 | 2.40e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 58 | 0.2069 | 0.2814 | 0.8 | 0.62 | 1.04 | 9.81e-02 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
