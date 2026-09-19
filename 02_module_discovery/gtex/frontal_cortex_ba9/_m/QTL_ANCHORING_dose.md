# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 7881 | 0.2223 | 0.2046 | 0.86 | 0.82 | 0.9 | 4.82e-11 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 5360 | 0.2155 | 0.2144 | 0.85 | 0.81 | 0.89 | 1.40e-12 | poisson_matched_standard |
| sQTL | go_invisible_modules | 702 | 0.235 | 0.2138 | 0.88 | 0.8 | 0.97 | 1.23e-02 | poisson_matched_standard |
| sQTL | go_visible_modules | 4658 | 0.2125 | 0.2161 | 0.86 | 0.82 | 0.91 | 1.82e-09 | poisson_matched_standard |
| eQTL | all_modules | 9201 | 0.5089 | 0.5525 | 0.87 | 0.83 | 0.91 | 2.85e-10 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 6299 | 0.5069 | 0.5427 | 0.87 | 0.83 | 0.91 | 1.50e-08 | poisson_matched_standard |
| eQTL | go_invisible_modules | 775 | 0.5097 | 0.5311 | 0.92 | 0.82 | 1.02 | 1.23e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 5524 | 0.5065 | 0.5407 | 0.88 | 0.84 | 0.92 | 2.43e-07 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
