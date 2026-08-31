# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3628 | 0.2412 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3013 | 0.2443 | 0.226 | 0.96 | 0.89 | 1.05 | 3.67e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1140 | 0.2807 | 0.2231 | 1.04 | 0.98 | 1.12 | 2.12e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 1873 | 0.2221 | 0.2615 | 0.94 | 0.89 | 1.0 | 6.68e-02 | ols_rankint_matched_standard |
| eQTL | all_modules | 4298 | 0.5209 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3488 | 0.531 | 0.4778 | 1.17 | 1.08 | 1.26 | 6.06e-05 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1298 | 0.5593 | 0.5043 | 1.12 | 1.05 | 1.2 | 5.74e-04 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2190 | 0.5142 | 0.528 | 1.0 | 0.94 | 1.06 | 9.96e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
