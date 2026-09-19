# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6105 | 0.1242 | 0.1191 | 0.97 | 0.93 | 1.0 | 5.44e-02 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 50 | 0.02 | 0.1218 | 0.82 | 0.63 | 1.08 | 1.68e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 50 | 0.02 | 0.1218 | 0.82 | 0.63 | 1.08 | 1.68e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 7044 | 0.2656 | 0.2998 | 0.9 | 0.88 | 0.93 | 6.90e-11 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 68 | 0.2059 | 0.2864 | 0.82 | 0.65 | 1.04 | 1.07e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 68 | 0.2059 | 0.2864 | 0.82 | 0.65 | 1.04 | 1.07e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
