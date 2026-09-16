# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3980 | 0.2317 | 0.2219 | 0.86 | 0.82 | 0.91 | 1.37e-09 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 2551 | 0.2148 | 0.2271 | 0.84 | 0.8 | 0.89 | 1.94e-09 | poisson_matched_standard |
| sQTL | go_invisible_modules | 511 | 0.2329 | 0.2245 | 0.97 | 0.87 | 1.08 | 6.17e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 2040 | 0.2103 | 0.2273 | 0.82 | 0.77 | 0.87 | 3.89e-10 | poisson_matched_standard |
| eQTL | all_modules | 4510 | 0.5182 | 0.5449 | 0.86 | 0.82 | 0.91 | 1.86e-08 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 2904 | 0.5017 | 0.5452 | 0.84 | 0.79 | 0.9 | 7.87e-08 | poisson_matched_standard |
| eQTL | go_invisible_modules | 568 | 0.5123 | 0.5391 | 0.88 | 0.77 | 1.0 | 5.48e-02 | poisson_matched_standard |
| eQTL | go_visible_modules | 2336 | 0.4991 | 0.544 | 0.84 | 0.78 | 0.9 | 1.15e-06 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
