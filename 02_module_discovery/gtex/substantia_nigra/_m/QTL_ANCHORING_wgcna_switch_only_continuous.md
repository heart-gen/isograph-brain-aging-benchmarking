# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5506 | 0.1095 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4825 | 0.1086 | 0.116 | 0.87 | 0.8 | 0.94 | 3.73e-04 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 148 | 0.0676 | 0.1107 | 0.96 | 0.82 | 1.13 | 6.60e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 4677 | 0.1099 | 0.1074 | 0.89 | 0.83 | 0.96 | 2.12e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 6721 | 0.257 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5778 | 0.2484 | 0.3097 | 0.88 | 0.83 | 0.95 | 4.50e-04 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 158 | 0.1962 | 0.2584 | 0.93 | 0.79 | 1.09 | 3.54e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 5620 | 0.2498 | 0.2934 | 0.91 | 0.85 | 0.97 | 3.59e-03 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
