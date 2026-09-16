# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4196 | 0.1118 | 0.126 | 0.73 | 0.65 | 0.82 | 1.78e-07 | logit_matched_standard |
| sQTL | pheno_sig_modules | 2745 | 0.0991 | 0.1274 | 0.68 | 0.59 | 0.79 | 1.50e-07 | logit_matched_standard |
| sQTL | go_invisible_modules | 631 | 0.1347 | 0.1208 | 0.9 | 0.71 | 1.15 | 4.04e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 2114 | 0.0885 | 0.1278 | 0.64 | 0.54 | 0.75 | 7.85e-08 | logit_matched_standard |
| eQTL | all_modules | 4823 | 0.2451 | 0.3016 | 0.73 | 0.67 | 0.79 | 4.65e-16 | logit_matched_standard |
| eQTL | pheno_sig_modules | 3149 | 0.2232 | 0.2998 | 0.67 | 0.61 | 0.73 | 2.85e-18 | logit_matched_standard |
| eQTL | go_invisible_modules | 697 | 0.2539 | 0.2874 | 0.79 | 0.66 | 0.94 | 8.45e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 2452 | 0.2145 | 0.2977 | 0.65 | 0.59 | 0.72 | 2.90e-16 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
