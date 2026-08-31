# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Spinal_cord_cervical_c-1 xQTL (gtex-aging/spinal_cord_cervical_c_1)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region spinal_cord_cervical_c_1`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3698 | 0.146 | 0.1574 | 0.77 | 0.72 | 0.82 | 2.63e-17 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 46 | 0.1304 | 0.1543 | 0.9 | 0.55 | 1.47 | 6.70e-01 | poisson_matched_standard |
| sQTL | go_invisible_modules | 18 | 0.1667 | 0.1542 | 0.87 | 0.39 | 1.94 | 7.36e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 28 | 0.1071 | 0.1543 | 0.92 | 0.49 | 1.7 | 7.82e-01 | poisson_matched_standard |
| eQTL | all_modules | 4761 | 0.3398 | 0.3789 | 0.79 | 0.74 | 0.84 | 3.29e-13 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 59 | 0.4407 | 0.3685 | 1.11 | 0.71 | 1.74 | 6.59e-01 | poisson_matched_standard |
| eQTL | go_invisible_modules | 21 | 0.5238 | 0.3686 | 0.96 | 0.43 | 2.15 | 9.29e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 38 | 0.3947 | 0.3687 | 1.19 | 0.69 | 2.05 | 5.36e-01 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
