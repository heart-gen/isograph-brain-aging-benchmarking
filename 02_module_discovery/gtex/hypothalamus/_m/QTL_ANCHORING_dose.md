# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4233 | 0.1857 | 0.1808 | 0.82 | 0.78 | 0.87 | 2.37e-12 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 482 | 0.1349 | 0.184 | 0.66 | 0.56 | 0.78 | 9.55e-07 | poisson_matched_standard |
| sQTL | go_invisible_modules | 216 | 0.1806 | 0.1823 | 0.92 | 0.75 | 1.12 | 4.13e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 266 | 0.0977 | 0.1839 | 0.41 | 0.3 | 0.55 | 3.41e-09 | poisson_matched_standard |
| eQTL | all_modules | 5183 | 0.3556 | 0.3985 | 0.8 | 0.76 | 0.85 | 9.40e-13 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 608 | 0.3026 | 0.3895 | 0.67 | 0.57 | 0.79 | 2.86e-06 | poisson_matched_standard |
| eQTL | go_invisible_modules | 234 | 0.3333 | 0.3874 | 0.88 | 0.7 | 1.11 | 2.88e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 374 | 0.2834 | 0.3888 | 0.55 | 0.44 | 0.69 | 4.76e-07 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
