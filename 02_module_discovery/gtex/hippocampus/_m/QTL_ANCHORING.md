# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6678 | 0.1628 | 0.1709 | 0.75 | 0.68 | 0.82 | 2.59e-09 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1165 | 0.1279 | 0.1705 | 0.7 | 0.58 | 0.84 | 1.44e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 40 | 0.225 | 0.1666 | 1.41 | 0.65 | 3.04 | 3.85e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1125 | 0.1244 | 0.1707 | 0.68 | 0.56 | 0.82 | 5.39e-05 | logit_matched_standard |
| eQTL | all_modules | 7780 | 0.3695 | 0.4059 | 0.82 | 0.77 | 0.87 | 2.01e-10 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1458 | 0.2867 | 0.3992 | 0.62 | 0.55 | 0.7 | 3.86e-15 | logit_matched_standard |
| eQTL | go_invisible_modules | 43 | 0.4419 | 0.3899 | 1.22 | 0.67 | 2.24 | 5.12e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1415 | 0.282 | 0.3993 | 0.61 | 0.54 | 0.69 | 6.36e-16 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
