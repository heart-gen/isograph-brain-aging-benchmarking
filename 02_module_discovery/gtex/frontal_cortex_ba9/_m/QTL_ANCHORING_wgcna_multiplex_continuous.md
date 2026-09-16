# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13359 | 0.2129 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 8234 | 0.2253 | 0.193 | 1.01 | 0.98 | 1.05 | 4.30e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 3089 | 0.2069 | 0.2147 | 0.94 | 0.91 | 0.98 | 3.80e-03 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 6649 | 0.2355 | 0.1905 | 1.06 | 1.02 | 1.09 | 6.72e-04 | ols_rankint_matched_standard |
| eQTL | all_modules | 17390 | 0.5278 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 10056 | 0.5369 | 0.5153 | 1.04 | 1.01 | 1.07 | 1.90e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 3463 | 0.5088 | 0.5325 | 0.92 | 0.88 | 0.95 | 5.69e-06 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 8213 | 0.5443 | 0.513 | 1.07 | 1.04 | 1.1 | 3.36e-06 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
