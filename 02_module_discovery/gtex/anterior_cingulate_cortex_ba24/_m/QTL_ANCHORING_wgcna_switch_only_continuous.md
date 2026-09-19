# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5996 | 0.1841 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4682 | 0.1773 | 0.2085 | 0.9 | 0.85 | 0.95 | 5.35e-04 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 4682 | 0.1773 | 0.2085 | 0.9 | 0.85 | 0.95 | 5.35e-04 | ols_rankint_matched_standard |
| eQTL | all_modules | 6573 | 0.3883 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5103 | 0.3837 | 0.4041 | 0.96 | 0.91 | 1.02 | 1.97e-01 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 5103 | 0.3837 | 0.4041 | 0.96 | 0.91 | 1.02 | 1.97e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
