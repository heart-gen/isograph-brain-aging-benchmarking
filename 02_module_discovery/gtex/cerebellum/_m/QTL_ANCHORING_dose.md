# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3643 | 0.3286 | 0.2649 | 0.92 | 0.89 | 0.96 | 3.42e-05 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 2101 | 0.3189 | 0.2757 | 0.87 | 0.83 | 0.91 | 2.45e-09 | poisson_matched_standard |
| sQTL | go_invisible_modules | 617 | 0.3679 | 0.2784 | 1.03 | 0.96 | 1.1 | 4.89e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 1484 | 0.2985 | 0.2806 | 0.81 | 0.76 | 0.85 | 6.64e-14 | poisson_matched_standard |
| eQTL | all_modules | 4058 | 0.6106 | 0.6131 | 0.96 | 0.92 | 1.01 | 8.33e-02 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 2346 | 0.5968 | 0.6149 | 0.93 | 0.87 | 0.98 | 1.10e-02 | poisson_matched_standard |
| eQTL | go_invisible_modules | 683 | 0.6047 | 0.6129 | 0.95 | 0.86 | 1.05 | 3.47e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 1663 | 0.5935 | 0.6145 | 0.92 | 0.86 | 0.99 | 1.96e-02 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
