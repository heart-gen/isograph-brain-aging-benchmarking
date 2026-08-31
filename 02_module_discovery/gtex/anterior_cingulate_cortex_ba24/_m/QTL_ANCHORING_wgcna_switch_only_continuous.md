# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6489 | 0.1581 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4412 | 0.1632 | 0.1473 | 1.0 | 0.95 | 1.06 | 8.77e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 95 | 0.2211 | 0.1572 | 1.29 | 1.05 | 1.57 | 1.35e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 4317 | 0.1619 | 0.1506 | 0.99 | 0.94 | 1.04 | 6.36e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 7671 | 0.3847 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5189 | 0.3849 | 0.3844 | 1.03 | 0.98 | 1.08 | 2.22e-01 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 112 | 0.3839 | 0.3847 | 1.07 | 0.89 | 1.29 | 4.62e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 5077 | 0.3849 | 0.3843 | 1.02 | 0.98 | 1.07 | 3.06e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
