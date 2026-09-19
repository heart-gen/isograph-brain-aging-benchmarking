# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5816 | 0.2395 | 0.2137 | 0.93 | 0.85 | 1.01 | 9.84e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2712 | 0.2297 | 0.2235 | 0.88 | 0.79 | 0.98 | 1.89e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 818 | 0.2726 | 0.2217 | 1.01 | 0.85 | 1.2 | 8.99e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1894 | 0.2112 | 0.227 | 0.84 | 0.74 | 0.95 | 4.93e-03 | logit_matched_standard |
| eQTL | all_modules | 6594 | 0.5253 | 0.5456 | 0.89 | 0.83 | 0.94 | 1.47e-04 | logit_matched_standard |
| eQTL | pheno_sig_modules | 3097 | 0.5082 | 0.5444 | 0.85 | 0.78 | 0.92 | 3.98e-05 | logit_matched_standard |
| eQTL | go_invisible_modules | 877 | 0.5564 | 0.5373 | 1.02 | 0.89 | 1.17 | 8.21e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 2220 | 0.4892 | 0.5451 | 0.8 | 0.73 | 0.88 | 1.28e-06 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
