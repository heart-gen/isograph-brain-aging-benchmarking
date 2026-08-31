# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3065 | 0.1171 | 0.1224 | 0.75 | 0.7 | 0.81 | 7.77e-13 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 1159 | 0.0984 | 0.1234 | 0.68 | 0.6 | 0.77 | 2.13e-09 | poisson_matched_standard |
| sQTL | go_invisible_modules | 515 | 0.0971 | 0.1222 | 0.71 | 0.59 | 0.86 | 4.24e-04 | poisson_matched_standard |
| sQTL | go_visible_modules | 644 | 0.0994 | 0.1223 | 0.69 | 0.59 | 0.81 | 5.33e-06 | poisson_matched_standard |
| eQTL | all_modules | 3840 | 0.2531 | 0.2949 | 0.78 | 0.72 | 0.84 | 4.89e-10 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 1507 | 0.2236 | 0.2917 | 0.71 | 0.62 | 0.8 | 8.60e-08 | poisson_matched_standard |
| eQTL | go_invisible_modules | 575 | 0.2243 | 0.2881 | 0.75 | 0.62 | 0.92 | 4.61e-03 | poisson_matched_standard |
| eQTL | go_visible_modules | 932 | 0.2232 | 0.2894 | 0.7 | 0.6 | 0.82 | 1.07e-05 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
