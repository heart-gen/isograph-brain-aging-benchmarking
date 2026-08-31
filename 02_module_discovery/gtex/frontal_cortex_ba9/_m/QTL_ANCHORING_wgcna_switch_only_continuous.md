# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6553 | 0.2109 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 5113 | 0.2132 | 0.2028 | 0.99 | 0.94 | 1.05 | 8.28e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 771 | 0.2296 | 0.2084 | 1.06 | 0.98 | 1.13 | 1.49e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 4342 | 0.2103 | 0.2121 | 0.97 | 0.92 | 1.02 | 2.40e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 7651 | 0.5109 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5877 | 0.5125 | 0.5056 | 1.05 | 1.0 | 1.11 | 6.55e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 858 | 0.4965 | 0.5127 | 0.99 | 0.92 | 1.06 | 8.19e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 5019 | 0.5152 | 0.5027 | 1.04 | 1.0 | 1.09 | 7.37e-02 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
