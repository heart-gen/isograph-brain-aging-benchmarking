# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 9878 | 0.2671 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7565 | 0.2619 | 0.284 | 0.93 | 0.89 | 0.97 | 6.55e-04 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 301 | 0.3621 | 0.2641 | 1.19 | 1.07 | 1.33 | 1.46e-03 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 7304 | 0.259 | 0.2898 | 0.91 | 0.87 | 0.95 | 1.09e-05 | ols_rankint_matched_standard |
| eQTL | all_modules | 12647 | 0.5774 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 9186 | 0.5837 | 0.5605 | 1.05 | 1.01 | 1.09 | 2.54e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 384 | 0.7318 | 0.5725 | 1.45 | 1.31 | 1.61 | 4.49e-13 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 8843 | 0.578 | 0.576 | 1.0 | 0.96 | 1.03 | 8.09e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
