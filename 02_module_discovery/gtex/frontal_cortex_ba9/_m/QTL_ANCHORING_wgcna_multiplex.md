# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13359 | 0.2129 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 8234 | 0.2253 | 0.193 | 1.08 | 0.99 | 1.18 | 9.15e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 3089 | 0.2069 | 0.2147 | 0.82 | 0.74 | 0.91 | 2.01e-04 | logit_matched_standard |
| sQTL | go_visible_modules | 6649 | 0.2355 | 0.1905 | 1.22 | 1.12 | 1.33 | 5.41e-06 | logit_matched_standard |
| eQTL | all_modules | 17390 | 0.5278 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 10056 | 0.5369 | 0.5153 | 1.09 | 1.03 | 1.16 | 5.79e-03 | logit_matched_standard |
| eQTL | go_invisible_modules | 3463 | 0.5088 | 0.5325 | 0.89 | 0.83 | 0.96 | 2.59e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 8213 | 0.5443 | 0.513 | 1.14 | 1.07 | 1.21 | 3.23e-05 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
