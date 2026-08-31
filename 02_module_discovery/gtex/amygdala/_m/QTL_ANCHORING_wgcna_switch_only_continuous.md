# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4858 | 0.1136 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3751 | 0.109 | 0.1292 | 0.96 | 0.9 | 1.03 | 2.51e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1982 | 0.1054 | 0.1193 | 0.97 | 0.92 | 1.02 | 2.66e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 1769 | 0.1131 | 0.114 | 1.0 | 0.95 | 1.06 | 8.89e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 6069 | 0.2595 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4591 | 0.2566 | 0.2686 | 0.94 | 0.89 | 1.0 | 5.19e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 2247 | 0.2555 | 0.2619 | 0.99 | 0.94 | 1.04 | 7.34e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2344 | 0.2577 | 0.2607 | 0.96 | 0.92 | 1.02 | 1.70e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
