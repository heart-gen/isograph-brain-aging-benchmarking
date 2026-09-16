# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4995 | 0.1562 | 0.1794 | 0.92 | 0.88 | 0.95 | 4.47e-07 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 2786 | 0.1439 | 0.1778 | 0.89 | 0.86 | 0.93 | 8.99e-08 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 425 | 0.2188 | 0.1692 | 1.02 | 0.93 | 1.12 | 6.97e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 2361 | 0.1305 | 0.1793 | 0.88 | 0.84 | 0.92 | 3.89e-09 | ols_rankint_matched_standard |
| eQTL | all_modules | 5823 | 0.3679 | 0.4249 | 0.87 | 0.85 | 0.9 | 2.09e-17 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 3354 | 0.3551 | 0.4181 | 0.86 | 0.83 | 0.9 | 2.28e-14 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 467 | 0.409 | 0.4061 | 1.0 | 0.91 | 1.1 | 9.99e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2887 | 0.3464 | 0.4178 | 0.85 | 0.81 | 0.88 | 4.69e-16 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
