# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6830 | 0.2151 | 0.2145 | 0.88 | 0.84 | 0.92 | 3.33e-08 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 5822 | 0.2133 | 0.2159 | 0.87 | 0.83 | 0.91 | 4.50e-09 | poisson_matched_standard |
| sQTL | go_invisible_modules | 2615 | 0.2528 | 0.2058 | 1.0 | 0.95 | 1.06 | 8.63e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 3207 | 0.1812 | 0.2251 | 0.81 | 0.76 | 0.86 | 5.12e-13 | poisson_matched_standard |
| eQTL | all_modules | 8211 | 0.5007 | 0.5542 | 0.78 | 0.74 | 0.81 | 1.50e-28 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 6983 | 0.5019 | 0.5475 | 0.8 | 0.76 | 0.84 | 3.08e-21 | poisson_matched_standard |
| eQTL | go_invisible_modules | 3138 | 0.5376 | 0.5285 | 1.0 | 0.94 | 1.06 | 9.74e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 3845 | 0.4728 | 0.5454 | 0.71 | 0.67 | 0.76 | 5.31e-29 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
