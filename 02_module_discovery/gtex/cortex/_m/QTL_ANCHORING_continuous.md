# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5816 | 0.2395 | 0.2137 | 0.96 | 0.93 | 0.99 | 1.05e-02 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 2712 | 0.2297 | 0.2235 | 0.93 | 0.89 | 0.97 | 3.34e-04 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 818 | 0.2726 | 0.2217 | 0.99 | 0.92 | 1.06 | 7.08e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 1894 | 0.2112 | 0.227 | 0.91 | 0.87 | 0.96 | 1.11e-04 | ols_rankint_matched_standard |
| eQTL | all_modules | 6594 | 0.5253 | 0.5456 | 0.92 | 0.9 | 0.95 | 2.99e-07 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 3097 | 0.5082 | 0.5444 | 0.9 | 0.87 | 0.94 | 1.50e-07 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 877 | 0.5564 | 0.5373 | 0.99 | 0.92 | 1.05 | 6.74e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2220 | 0.4892 | 0.5451 | 0.88 | 0.84 | 0.92 | 9.83e-09 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
