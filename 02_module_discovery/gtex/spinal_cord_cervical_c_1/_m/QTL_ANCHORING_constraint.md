# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Spinal_cord_cervical_c-1 xQTL (gtex-aging/spinal_cord_cervical_c_1)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region spinal_cord_cervical_c_1`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5026 | 0.1492 | 0.1305 | 0.96 | 0.86 | 1.07 | 4.66e-01 | logit_matched_standard |
| sQTL | all_modules | 5026 | 0.1492 | 0.1305 | 0.92 | 0.82 | 1.03 | 1.65e-01 | logit_matched_constraint |
| eQTL | all_modules | 5525 | 0.352 | 0.3386 | 1.0 | 0.93 | 1.08 | 9.30e-01 | logit_matched_standard |
| eQTL | all_modules | 5525 | 0.352 | 0.3386 | 0.99 | 0.91 | 1.07 | 7.35e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
