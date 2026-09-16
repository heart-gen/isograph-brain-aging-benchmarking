# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4607 | 0.153 | 0.1741 | 0.72 | 0.65 | 0.8 | 4.63e-10 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1406 | 0.1408 | 0.1699 | 0.74 | 0.63 | 0.87 | 3.37e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 33 | 0.0909 | 0.167 | 0.59 | 0.18 | 1.96 | 3.87e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1373 | 0.142 | 0.1696 | 0.75 | 0.63 | 0.88 | 4.92e-04 | logit_matched_standard |
| eQTL | all_modules | 5370 | 0.3512 | 0.4068 | 0.77 | 0.72 | 0.82 | 1.36e-14 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1731 | 0.3056 | 0.3991 | 0.67 | 0.6 | 0.75 | 3.10e-13 | logit_matched_standard |
| eQTL | go_invisible_modules | 40 | 0.3 | 0.3902 | 0.63 | 0.32 | 1.24 | 1.81e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1691 | 0.3057 | 0.3989 | 0.67 | 0.6 | 0.75 | 8.67e-13 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
