# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6995 | 0.1654 | 0.1766 | 0.94 | 0.91 | 0.97 | 1.04e-04 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 2628 | 0.1419 | 0.1778 | 0.89 | 0.85 | 0.93 | 2.66e-08 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 2628 | 0.1419 | 0.1778 | 0.89 | 0.85 | 0.93 | 2.66e-08 | ols_rankint_matched_standard |
| eQTL | all_modules | 8148 | 0.3798 | 0.4285 | 0.88 | 0.85 | 0.9 | 7.08e-18 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 3193 | 0.3539 | 0.4176 | 0.86 | 0.83 | 0.9 | 9.88e-14 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 3193 | 0.3539 | 0.4176 | 0.86 | 0.83 | 0.9 | 9.88e-14 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
