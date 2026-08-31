# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13376 | 0.2231 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7243 | 0.2369 | 0.2068 | 1.01 | 0.98 | 1.05 | 4.39e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1431 | 0.2739 | 0.217 | 1.03 | 0.98 | 1.09 | 2.68e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 6524 | 0.2333 | 0.2134 | 1.01 | 0.98 | 1.05 | 4.45e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 17847 | 0.5335 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 9134 | 0.5504 | 0.5159 | 1.09 | 1.06 | 1.12 | 4.67e-09 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1724 | 0.5505 | 0.5317 | 1.02 | 0.97 | 1.07 | 4.19e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 8236 | 0.5501 | 0.5193 | 1.08 | 1.05 | 1.11 | 1.17e-07 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
