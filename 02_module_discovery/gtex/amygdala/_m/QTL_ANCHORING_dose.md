# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6105 | 0.1242 | 0.1191 | 0.79 | 0.74 | 0.85 | 3.97e-12 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 50 | 0.02 | 0.1218 | 0.09 | 0.01 | 0.62 | 1.51e-02 | poisson_matched_standard |
| sQTL | go_visible_modules | 50 | 0.02 | 0.1218 | 0.09 | 0.01 | 0.62 | 1.51e-02 | poisson_matched_standard |
| eQTL | all_modules | 7044 | 0.2656 | 0.2998 | 0.8 | 0.75 | 0.85 | 1.64e-11 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 68 | 0.2059 | 0.2864 | 0.52 | 0.27 | 1.0 | 5.14e-02 | poisson_matched_standard |
| eQTL | go_visible_modules | 68 | 0.2059 | 0.2864 | 0.52 | 0.27 | 1.0 | 5.14e-02 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
