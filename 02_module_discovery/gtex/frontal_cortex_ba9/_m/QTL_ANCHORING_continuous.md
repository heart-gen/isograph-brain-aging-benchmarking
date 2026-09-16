# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5845 | 0.2087 | 0.2195 | 0.93 | 0.9 | 0.96 | 2.42e-05 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 4257 | 0.1976 | 0.2227 | 0.93 | 0.89 | 0.96 | 1.87e-05 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1801 | 0.2099 | 0.2156 | 0.94 | 0.9 | 0.99 | 1.19e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 2456 | 0.1885 | 0.2207 | 0.94 | 0.9 | 0.98 | 3.52e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 6846 | 0.4963 | 0.551 | 0.88 | 0.85 | 0.9 | 1.51e-17 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 5056 | 0.4941 | 0.5443 | 0.89 | 0.86 | 0.92 | 5.12e-12 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1985 | 0.4922 | 0.5349 | 0.91 | 0.86 | 0.95 | 3.28e-05 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 3071 | 0.4953 | 0.5374 | 0.91 | 0.87 | 0.95 | 1.79e-06 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
