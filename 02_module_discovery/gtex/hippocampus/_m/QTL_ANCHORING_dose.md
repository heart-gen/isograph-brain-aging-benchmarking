# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4607 | 0.153 | 0.1741 | 0.73 | 0.69 | 0.78 | 3.28e-26 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 1406 | 0.1408 | 0.1699 | 0.77 | 0.7 | 0.84 | 1.77e-08 | poisson_matched_standard |
| sQTL | go_invisible_modules | 33 | 0.0909 | 0.167 | 0.83 | 0.42 | 1.66 | 6.03e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 1373 | 0.142 | 0.1696 | 0.77 | 0.7 | 0.84 | 2.10e-08 | poisson_matched_standard |
| eQTL | all_modules | 5370 | 0.3512 | 0.4068 | 0.78 | 0.74 | 0.83 | 6.37e-16 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 1731 | 0.3056 | 0.3991 | 0.69 | 0.63 | 0.77 | 1.26e-12 | poisson_matched_standard |
| eQTL | go_invisible_modules | 40 | 0.3 | 0.3902 | 0.5 | 0.24 | 1.05 | 6.76e-02 | poisson_matched_standard |
| eQTL | go_visible_modules | 1691 | 0.3057 | 0.3989 | 0.7 | 0.63 | 0.77 | 6.40e-12 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
