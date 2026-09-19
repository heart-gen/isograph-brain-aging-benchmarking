# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5175 | 0.3277 | 0.2531 | 1.01 | 0.93 | 1.09 | 8.79e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2869 | 0.3283 | 0.2698 | 0.98 | 0.89 | 1.08 | 7.00e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1001 | 0.3666 | 0.2757 | 1.07 | 0.93 | 1.24 | 3.59e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1868 | 0.3078 | 0.2785 | 0.93 | 0.83 | 1.05 | 2.39e-01 | logit_matched_standard |
| eQTL | all_modules | 5785 | 0.6145 | 0.6117 | 0.97 | 0.9 | 1.03 | 2.94e-01 | logit_matched_standard |
| eQTL | pheno_sig_modules | 3187 | 0.6056 | 0.614 | 0.93 | 0.86 | 1.01 | 6.68e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 1098 | 0.6138 | 0.6125 | 0.98 | 0.86 | 1.11 | 7.00e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 2089 | 0.6012 | 0.614 | 0.91 | 0.83 | 1.0 | 5.83e-02 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
