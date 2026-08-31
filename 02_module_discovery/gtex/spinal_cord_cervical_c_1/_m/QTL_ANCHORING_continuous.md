# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Spinal_cord_cervical_c-1 xQTL (gtex-aging/spinal_cord_cervical_c_1)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region spinal_cord_cervical_c_1`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3698 | 0.146 | 0.1574 | 0.94 | 0.9 | 0.98 | 1.21e-03 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 46 | 0.1304 | 0.1543 | 0.89 | 0.67 | 1.19 | 4.39e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 18 | 0.1667 | 0.1542 | 1.22 | 0.78 | 1.93 | 3.80e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 28 | 0.1071 | 0.1543 | 0.73 | 0.51 | 1.05 | 9.01e-02 | ols_rankint_matched_standard |
| eQTL | all_modules | 4761 | 0.3398 | 0.3789 | 0.88 | 0.85 | 0.91 | 1.36e-14 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 59 | 0.4407 | 0.3685 | 1.09 | 0.84 | 1.4 | 5.19e-01 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 21 | 0.5238 | 0.3686 | 1.21 | 0.79 | 1.85 | 3.87e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 38 | 0.3947 | 0.3687 | 1.03 | 0.75 | 1.41 | 8.73e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
