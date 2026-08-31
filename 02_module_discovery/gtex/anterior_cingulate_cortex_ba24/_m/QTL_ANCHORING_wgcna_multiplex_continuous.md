# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13291 | 0.1697 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 10423 | 0.1659 | 0.1834 | 0.94 | 0.9 | 0.98 | 2.74e-03 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1744 | 0.1927 | 0.1662 | 0.99 | 0.94 | 1.04 | 6.69e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 10086 | 0.1648 | 0.185 | 0.94 | 0.91 | 0.98 | 2.70e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 17572 | 0.4033 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 13566 | 0.4017 | 0.4086 | 1.01 | 0.97 | 1.04 | 7.79e-01 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 2136 | 0.4167 | 0.4015 | 1.0 | 0.95 | 1.04 | 8.68e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 12991 | 0.3975 | 0.4198 | 0.98 | 0.95 | 1.02 | 3.23e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
