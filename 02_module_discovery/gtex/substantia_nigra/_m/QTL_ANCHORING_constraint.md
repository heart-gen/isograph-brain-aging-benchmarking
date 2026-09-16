# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3868 | 0.1122 | 0.1096 | 0.9 | 0.79 | 1.02 | 9.72e-02 | logit_matched_standard |
| sQTL | all_modules | 3868 | 0.1122 | 0.1096 | 0.83 | 0.73 | 0.95 | 6.16e-03 | logit_matched_constraint |
| eQTL | all_modules | 4240 | 0.2356 | 0.2573 | 0.86 | 0.79 | 0.94 | 6.96e-04 | logit_matched_standard |
| eQTL | all_modules | 4240 | 0.2356 | 0.2573 | 0.88 | 0.8 | 0.96 | 4.20e-03 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
