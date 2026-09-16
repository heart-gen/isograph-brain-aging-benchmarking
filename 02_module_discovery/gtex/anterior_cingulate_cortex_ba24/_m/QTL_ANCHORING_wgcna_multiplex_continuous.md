# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13236 | 0.1698 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 10490 | 0.165 | 0.1883 | 0.93 | 0.89 | 0.97 | 7.53e-04 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 1754 | 0.1973 | 0.1657 | 1.03 | 0.98 | 1.08 | 2.96e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 10175 | 0.1631 | 0.1921 | 0.92 | 0.88 | 0.96 | 2.96e-05 | ols_rankint_matched_standard |
| eQTL | all_modules | 17280 | 0.4031 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 13468 | 0.3956 | 0.4294 | 0.96 | 0.92 | 0.99 | 1.42e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1927 | 0.4001 | 0.4034 | 0.95 | 0.91 | 1.0 | 3.99e-02 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 13095 | 0.3939 | 0.4318 | 0.95 | 0.92 | 0.99 | 5.38e-03 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
