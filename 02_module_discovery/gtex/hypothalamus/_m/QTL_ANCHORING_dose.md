# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6581 | 0.1854 | 0.1796 | 0.79 | 0.75 | 0.83 | 1.09e-20 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 1111 | 0.18 | 0.1826 | 0.83 | 0.76 | 0.91 | 9.19e-05 | poisson_matched_standard |
| sQTL | go_invisible_modules | 1039 | 0.1829 | 0.1823 | 0.81 | 0.74 | 0.9 | 2.25e-05 | poisson_matched_standard |
| sQTL | go_visible_modules | 72 | 0.1389 | 0.1826 | 1.15 | 0.84 | 1.57 | 3.92e-01 | poisson_matched_standard |
| eQTL | all_modules | 7457 | 0.3731 | 0.3976 | 0.87 | 0.83 | 0.92 | 3.22e-07 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 1258 | 0.3498 | 0.3904 | 0.92 | 0.83 | 1.02 | 1.15e-01 | poisson_matched_standard |
| eQTL | go_invisible_modules | 1159 | 0.3529 | 0.39 | 0.93 | 0.84 | 1.03 | 1.79e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 99 | 0.3131 | 0.3881 | 0.82 | 0.56 | 1.21 | 3.17e-01 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
