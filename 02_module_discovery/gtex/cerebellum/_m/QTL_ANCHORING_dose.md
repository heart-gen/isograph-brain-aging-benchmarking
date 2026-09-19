# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5175 | 0.3277 | 0.2531 | 0.91 | 0.88 | 0.95 | 1.82e-06 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 2869 | 0.3283 | 0.2698 | 0.92 | 0.88 | 0.96 | 5.14e-05 | poisson_matched_standard |
| sQTL | go_invisible_modules | 1001 | 0.3666 | 0.2757 | 1.03 | 0.97 | 1.09 | 3.31e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 1868 | 0.3078 | 0.2785 | 0.87 | 0.82 | 0.91 | 2.51e-08 | poisson_matched_standard |
| eQTL | all_modules | 5785 | 0.6145 | 0.6117 | 0.96 | 0.92 | 1.0 | 6.92e-02 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 3187 | 0.6056 | 0.614 | 0.95 | 0.9 | 1.0 | 3.75e-02 | poisson_matched_standard |
| eQTL | go_invisible_modules | 1098 | 0.6138 | 0.6125 | 0.98 | 0.9 | 1.06 | 5.91e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 2089 | 0.6012 | 0.614 | 0.94 | 0.88 | 1.0 | 3.81e-02 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
