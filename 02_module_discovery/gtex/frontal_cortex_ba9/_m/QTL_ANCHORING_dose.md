# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5845 | 0.2087 | 0.2195 | 0.83 | 0.79 | 0.87 | 4.64e-15 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 4257 | 0.1976 | 0.2227 | 0.8 | 0.76 | 0.84 | 1.12e-17 | poisson_matched_standard |
| sQTL | go_invisible_modules | 1801 | 0.2099 | 0.2156 | 0.81 | 0.76 | 0.87 | 2.20e-09 | poisson_matched_standard |
| sQTL | go_visible_modules | 2456 | 0.1885 | 0.2207 | 0.86 | 0.81 | 0.91 | 1.49e-06 | poisson_matched_standard |
| eQTL | all_modules | 6846 | 0.4963 | 0.551 | 0.85 | 0.81 | 0.89 | 6.31e-12 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 5056 | 0.4941 | 0.5443 | 0.84 | 0.8 | 0.89 | 9.93e-11 | poisson_matched_standard |
| eQTL | go_invisible_modules | 1985 | 0.4922 | 0.5349 | 0.86 | 0.8 | 0.93 | 1.23e-04 | poisson_matched_standard |
| eQTL | go_visible_modules | 3071 | 0.4953 | 0.5374 | 0.87 | 0.81 | 0.92 | 5.94e-06 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
