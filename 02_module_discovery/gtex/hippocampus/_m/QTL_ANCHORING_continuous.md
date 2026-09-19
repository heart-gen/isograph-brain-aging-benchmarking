# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6678 | 0.1628 | 0.1709 | 0.92 | 0.89 | 0.95 | 7.73e-07 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 1165 | 0.1279 | 0.1705 | 0.89 | 0.84 | 0.94 | 5.85e-05 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 40 | 0.225 | 0.1666 | 0.98 | 0.72 | 1.33 | 8.94e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 1125 | 0.1244 | 0.1707 | 0.88 | 0.83 | 0.94 | 4.98e-05 | ols_rankint_matched_standard |
| eQTL | all_modules | 7780 | 0.3695 | 0.4059 | 0.91 | 0.88 | 0.94 | 2.81e-10 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 1458 | 0.2867 | 0.3992 | 0.82 | 0.77 | 0.86 | 7.73e-14 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 43 | 0.4419 | 0.3899 | 0.86 | 0.64 | 1.16 | 3.23e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 1415 | 0.282 | 0.3993 | 0.82 | 0.77 | 0.86 | 1.38e-13 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
