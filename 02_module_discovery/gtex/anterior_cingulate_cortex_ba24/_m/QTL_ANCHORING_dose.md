# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6995 | 0.1654 | 0.1766 | 0.78 | 0.74 | 0.82 | 2.75e-21 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 2628 | 0.1419 | 0.1778 | 0.75 | 0.7 | 0.81 | 4.61e-15 | poisson_matched_standard |
| sQTL | go_visible_modules | 2628 | 0.1419 | 0.1778 | 0.75 | 0.7 | 0.81 | 4.61e-15 | poisson_matched_standard |
| eQTL | all_modules | 8148 | 0.3798 | 0.4285 | 0.79 | 0.75 | 0.83 | 4.08e-18 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 3193 | 0.3539 | 0.4176 | 0.75 | 0.69 | 0.8 | 5.09e-14 | poisson_matched_standard |
| eQTL | go_visible_modules | 3193 | 0.3539 | 0.4176 | 0.75 | 0.69 | 0.8 | 5.09e-14 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
