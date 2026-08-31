# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2527 | 0.2489 | 0.2196 | 0.95 | 0.91 | 1.01 | 8.98e-02 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 1380 | 0.2312 | 0.2243 | 0.82 | 0.76 | 0.88 | 4.11e-08 | poisson_matched_standard |
| sQTL | go_invisible_modules | 824 | 0.2536 | 0.2232 | 0.77 | 0.7 | 0.84 | 6.94e-09 | poisson_matched_standard |
| sQTL | go_visible_modules | 556 | 0.1978 | 0.2262 | 0.94 | 0.84 | 1.05 | 2.91e-01 | poisson_matched_standard |
| eQTL | all_modules | 3164 | 0.5035 | 0.5442 | 0.83 | 0.78 | 0.88 | 1.77e-09 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 1676 | 0.4875 | 0.5422 | 0.78 | 0.72 | 0.85 | 4.19e-09 | poisson_matched_standard |
| eQTL | go_invisible_modules | 879 | 0.5006 | 0.539 | 0.82 | 0.73 | 0.91 | 2.47e-04 | poisson_matched_standard |
| eQTL | go_visible_modules | 797 | 0.473 | 0.5401 | 0.76 | 0.68 | 0.86 | 8.66e-06 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
