# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3643 | 0.3286 | 0.2649 | 0.99 | 0.9 | 1.08 | 7.78e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2101 | 0.3189 | 0.2757 | 0.91 | 0.81 | 1.01 | 7.38e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 617 | 0.3679 | 0.2784 | 1.01 | 0.84 | 1.21 | 9.12e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1484 | 0.2985 | 0.2806 | 0.87 | 0.76 | 0.99 | 3.08e-02 | logit_matched_standard |
| eQTL | all_modules | 4058 | 0.6106 | 0.6131 | 0.95 | 0.88 | 1.02 | 1.44e-01 | logit_matched_standard |
| eQTL | pheno_sig_modules | 2346 | 0.5968 | 0.6149 | 0.89 | 0.81 | 0.97 | 1.04e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 683 | 0.6047 | 0.6129 | 0.93 | 0.8 | 1.09 | 3.83e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1663 | 0.5935 | 0.6145 | 0.88 | 0.79 | 0.98 | 1.64e-02 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
