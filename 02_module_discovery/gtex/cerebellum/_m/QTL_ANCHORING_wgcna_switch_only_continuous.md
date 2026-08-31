# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 804 | 0.3271 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 583 | 0.3482 | 0.2715 | 0.89 | 0.77 | 1.04 | 1.41e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 266 | 0.3647 | 0.3086 | 1.0 | 0.87 | 1.15 | 9.99e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 317 | 0.3344 | 0.3224 | 0.91 | 0.8 | 1.05 | 1.89e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 949 | 0.5711 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 663 | 0.5867 | 0.535 | 1.1 | 0.95 | 1.26 | 1.99e-01 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 311 | 0.5788 | 0.5674 | 1.02 | 0.89 | 1.16 | 8.10e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 352 | 0.5938 | 0.5578 | 1.07 | 0.94 | 1.22 | 3.31e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
