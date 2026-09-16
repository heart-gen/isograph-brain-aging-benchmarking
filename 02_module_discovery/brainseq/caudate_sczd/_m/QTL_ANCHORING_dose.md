# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4607 | 0.2129 | 0.2306 | 0.86 | 0.82 | 0.9 | 1.97e-09 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 417 | 0.2518 | 0.2238 | 0.96 | 0.85 | 1.09 | 5.46e-01 | poisson_matched_standard |
| sQTL | go_invisible_modules | 380 | 0.2553 | 0.2237 | 0.93 | 0.82 | 1.07 | 3.15e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 37 | 0.2162 | 0.2246 | 1.34 | 0.89 | 2.01 | 1.66e-01 | poisson_matched_standard |
| eQTL | all_modules | 5816 | 0.4713 | 0.5494 | 0.74 | 0.7 | 0.78 | 2.64e-30 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 458 | 0.548 | 0.5244 | 1.04 | 0.91 | 1.2 | 5.48e-01 | poisson_matched_standard |
| eQTL | go_invisible_modules | 415 | 0.5494 | 0.5244 | 1.05 | 0.91 | 1.22 | 4.89e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 43 | 0.5349 | 0.525 | 0.95 | 0.59 | 1.54 | 8.47e-01 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
