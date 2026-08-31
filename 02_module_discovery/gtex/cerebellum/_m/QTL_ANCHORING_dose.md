# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2355 | 0.3231 | 0.2735 | 0.86 | 0.82 | 0.9 | 8.84e-11 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 700 | 0.3071 | 0.281 | 0.87 | 0.81 | 0.94 | 2.71e-04 | poisson_matched_standard |
| sQTL | go_invisible_modules | 504 | 0.369 | 0.279 | 0.9 | 0.83 | 0.98 | 9.82e-03 | poisson_matched_standard |
| sQTL | go_visible_modules | 196 | 0.148 | 0.2845 | 0.71 | 0.57 | 0.88 | 1.64e-03 | poisson_matched_standard |
| eQTL | all_modules | 2782 | 0.5945 | 0.6133 | 0.86 | 0.81 | 0.91 | 1.08e-07 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 921 | 0.5765 | 0.6123 | 0.79 | 0.72 | 0.87 | 1.95e-06 | poisson_matched_standard |
| eQTL | go_invisible_modules | 552 | 0.6739 | 0.6086 | 1.02 | 0.92 | 1.14 | 6.97e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 369 | 0.4309 | 0.6141 | 0.46 | 0.38 | 0.57 | 1.47e-14 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
