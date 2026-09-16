# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4502 | 0.1766 | 0.1851 | 0.92 | 0.88 | 0.95 | 9.50e-07 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 338 | 0.1746 | 0.1825 | 0.92 | 0.83 | 1.03 | 1.41e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 338 | 0.1746 | 0.1825 | 0.92 | 0.83 | 1.03 | 1.41e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 5055 | 0.3589 | 0.3986 | 0.92 | 0.89 | 0.95 | 1.35e-07 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 373 | 0.3378 | 0.3887 | 0.89 | 0.81 | 0.99 | 3.10e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 373 | 0.3378 | 0.3887 | 0.89 | 0.81 | 0.99 | 3.10e-02 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
