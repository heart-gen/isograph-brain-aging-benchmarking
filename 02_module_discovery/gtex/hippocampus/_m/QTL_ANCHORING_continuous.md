# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5335 | 0.1693 | 0.1655 | 0.94 | 0.91 | 0.97 | 5.94e-04 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 2298 | 0.1658 | 0.1672 | 0.96 | 0.92 | 1.01 | 9.31e-02 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 945 | 0.1937 | 0.165 | 1.03 | 0.97 | 1.1 | 3.21e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 1353 | 0.1463 | 0.1693 | 0.92 | 0.87 | 0.97 | 3.20e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 6742 | 0.3558 | 0.4089 | 0.88 | 0.85 | 0.9 | 3.75e-17 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 3028 | 0.3322 | 0.4007 | 0.85 | 0.82 | 0.89 | 8.38e-16 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1094 | 0.372 | 0.3906 | 0.97 | 0.91 | 1.03 | 2.81e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 1934 | 0.3097 | 0.3988 | 0.81 | 0.77 | 0.85 | 7.13e-19 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
