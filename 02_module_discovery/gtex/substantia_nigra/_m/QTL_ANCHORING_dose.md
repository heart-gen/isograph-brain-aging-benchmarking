# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3886 | 0.1225 | 0.1252 | 0.74 | 0.69 | 0.79 | 4.17e-17 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 58 | 0.0862 | 0.1246 | 0.67 | 0.4 | 1.12 | 1.25e-01 | poisson_matched_standard |
| sQTL | go_invisible_modules | 30 | 0.1 | 0.1245 | 0.96 | 0.55 | 1.7 | 9.01e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 28 | 0.0714 | 0.1245 | 0.3 | 0.1 | 0.94 | 3.95e-02 | poisson_matched_standard |
| eQTL | all_modules | 4751 | 0.2572 | 0.2896 | 0.79 | 0.74 | 0.85 | 8.33e-10 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 98 | 0.1633 | 0.2818 | 0.52 | 0.29 | 0.91 | 2.24e-02 | poisson_matched_standard |
| eQTL | go_invisible_modules | 40 | 0.1 | 0.2815 | 0.42 | 0.16 | 1.12 | 8.44e-02 | poisson_matched_standard |
| eQTL | go_visible_modules | 58 | 0.2069 | 0.2814 | 0.58 | 0.29 | 1.17 | 1.29e-01 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
