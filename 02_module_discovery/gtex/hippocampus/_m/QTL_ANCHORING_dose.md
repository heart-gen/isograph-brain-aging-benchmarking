# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6678 | 0.1628 | 0.1709 | 0.73 | 0.69 | 0.77 | 5.48e-30 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 1165 | 0.1279 | 0.1705 | 0.78 | 0.7 | 0.86 | 2.24e-06 | poisson_matched_standard |
| sQTL | go_invisible_modules | 40 | 0.225 | 0.1666 | 1.43 | 0.96 | 2.12 | 7.56e-02 | poisson_matched_standard |
| sQTL | go_visible_modules | 1125 | 0.1244 | 0.1707 | 0.76 | 0.68 | 0.84 | 2.45e-07 | poisson_matched_standard |
| eQTL | all_modules | 7780 | 0.3695 | 0.4059 | 0.83 | 0.79 | 0.88 | 8.17e-12 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 1458 | 0.2867 | 0.3992 | 0.65 | 0.58 | 0.73 | 1.25e-13 | poisson_matched_standard |
| eQTL | go_invisible_modules | 43 | 0.4419 | 0.3899 | 1.18 | 0.73 | 1.9 | 4.90e-01 | poisson_matched_standard |
| eQTL | go_visible_modules | 1415 | 0.282 | 0.3993 | 0.63 | 0.56 | 0.71 | 2.21e-14 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
