# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3980 | 0.2317 | 0.2219 | 0.89 | 0.81 | 0.97 | 1.26e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2551 | 0.2148 | 0.2271 | 0.82 | 0.73 | 0.91 | 3.24e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 511 | 0.2329 | 0.2245 | 0.91 | 0.73 | 1.14 | 4.11e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 2040 | 0.2103 | 0.2273 | 0.81 | 0.71 | 0.91 | 5.11e-04 | logit_matched_standard |
| eQTL | all_modules | 4510 | 0.5182 | 0.5449 | 0.88 | 0.82 | 0.94 | 1.30e-04 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2904 | 0.5017 | 0.5452 | 0.83 | 0.77 | 0.9 | 6.94e-06 | logit_matched_standard |
| eQTL | go_invisible_modules | 568 | 0.5123 | 0.5391 | 0.86 | 0.73 | 1.02 | 8.15e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 2336 | 0.4991 | 0.544 | 0.84 | 0.77 | 0.91 | 6.36e-05 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
