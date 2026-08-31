# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13415 | 0.2141 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 8025 | 0.2204 | 0.2046 | 0.98 | 0.95 | 1.01 | 2.68e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1683 | 0.2436 | 0.2099 | 1.0 | 0.96 | 1.06 | 8.54e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 7320 | 0.2178 | 0.2097 | 0.99 | 0.96 | 1.02 | 4.98e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 17569 | 0.5273 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 10046 | 0.533 | 0.5197 | 1.03 | 1.0 | 1.06 | 5.40e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1883 | 0.5236 | 0.5278 | 0.96 | 0.91 | 1.0 | 6.03e-02 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 9249 | 0.5355 | 0.5183 | 1.04 | 1.01 | 1.07 | 4.69e-03 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
