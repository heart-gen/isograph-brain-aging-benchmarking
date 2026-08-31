# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4757 | 0.1682 | 0.172 | 0.87 | 0.82 | 0.92 | 9.69e-07 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 3236 | 0.1582 | 0.1746 | 0.82 | 0.77 | 0.87 | 1.11e-09 | poisson_matched_standard |
| sQTL | go_invisible_modules | 1652 | 0.181 | 0.1692 | 0.94 | 0.87 | 1.01 | 9.15e-02 | poisson_matched_standard |
| sQTL | go_visible_modules | 1584 | 0.1345 | 0.1755 | 0.74 | 0.67 | 0.81 | 1.33e-10 | poisson_matched_standard |
| eQTL | all_modules | 5799 | 0.3656 | 0.426 | 0.75 | 0.7 | 0.79 | 4.67e-22 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 3954 | 0.3488 | 0.4228 | 0.73 | 0.68 | 0.78 | 1.13e-18 | poisson_matched_standard |
| eQTL | go_invisible_modules | 1910 | 0.378 | 0.41 | 0.85 | 0.77 | 0.93 | 2.66e-04 | poisson_matched_standard |
| eQTL | go_visible_modules | 2044 | 0.3214 | 0.4175 | 0.68 | 0.61 | 0.75 | 2.43e-15 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
