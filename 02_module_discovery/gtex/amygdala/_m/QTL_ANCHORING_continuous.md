# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3065 | 0.1171 | 0.1224 | 0.94 | 0.91 | 0.98 | 4.89e-03 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 1159 | 0.0984 | 0.1234 | 0.93 | 0.88 | 0.99 | 2.57e-02 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 515 | 0.0971 | 0.1222 | 0.93 | 0.85 | 1.01 | 9.32e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 644 | 0.0994 | 0.1223 | 0.95 | 0.87 | 1.02 | 1.65e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 3840 | 0.2531 | 0.2949 | 0.89 | 0.86 | 0.92 | 8.67e-11 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 1507 | 0.2236 | 0.2917 | 0.84 | 0.8 | 0.89 | 1.85e-10 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 575 | 0.2243 | 0.2881 | 0.88 | 0.81 | 0.96 | 3.09e-03 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 932 | 0.2232 | 0.2894 | 0.83 | 0.78 | 0.89 | 2.57e-08 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
