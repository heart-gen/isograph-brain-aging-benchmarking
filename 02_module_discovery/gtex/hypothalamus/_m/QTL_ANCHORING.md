# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6581 | 0.1854 | 0.1796 | 0.85 | 0.77 | 0.93 | 3.50e-04 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1111 | 0.18 | 0.1826 | 0.84 | 0.71 | 1.0 | 4.50e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 1039 | 0.1829 | 0.1823 | 0.85 | 0.72 | 1.01 | 6.94e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 72 | 0.1389 | 0.1826 | 0.72 | 0.36 | 1.43 | 3.47e-01 | logit_matched_standard |
| eQTL | all_modules | 7457 | 0.3731 | 0.3976 | 0.87 | 0.81 | 0.92 | 6.18e-06 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1258 | 0.3498 | 0.3904 | 0.8 | 0.71 | 0.91 | 4.41e-04 | logit_matched_standard |
| eQTL | go_invisible_modules | 1159 | 0.3529 | 0.39 | 0.81 | 0.72 | 0.92 | 1.14e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 99 | 0.3131 | 0.3881 | 0.75 | 0.49 | 1.15 | 1.85e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
