# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13051 | 0.1211 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 5532 | 0.1215 | 0.1209 | 1.0 | 0.97 | 1.04 | 8.15e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 827 | 0.1197 | 0.1212 | 0.99 | 0.92 | 1.06 | 8.07e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 5142 | 0.1208 | 0.1214 | 1.0 | 0.97 | 1.04 | 9.88e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 17757 | 0.283 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7403 | 0.2645 | 0.2963 | 0.92 | 0.9 | 0.95 | 2.19e-07 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 989 | 0.2538 | 0.2848 | 0.94 | 0.88 | 1.0 | 4.43e-02 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 6933 | 0.2628 | 0.296 | 0.92 | 0.89 | 0.95 | 3.48e-08 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
