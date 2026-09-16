# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Cerebellar_Hemisphere xQTL (gtex-aging/cerebellar_hemisphere)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellar_hemisphere`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4456 | 0.309 | 0.2604 | 1.01 | 0.92 | 1.11 | 8.23e-01 | logit_matched_standard |
| sQTL | all_modules | 4456 | 0.309 | 0.2604 | 0.99 | 0.9 | 1.08 | 7.72e-01 | logit_matched_constraint |
| eQTL | all_modules | 4831 | 0.5957 | 0.5992 | 0.96 | 0.9 | 1.04 | 3.42e-01 | logit_matched_standard |
| eQTL | all_modules | 4831 | 0.5957 | 0.5992 | 0.91 | 0.85 | 0.99 | 2.21e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
