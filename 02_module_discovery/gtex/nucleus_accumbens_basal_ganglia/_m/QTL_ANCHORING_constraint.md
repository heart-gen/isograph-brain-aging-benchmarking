# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Nucleus_accumbens_basal_ganglia xQTL (gtex-aging/nucleus_accumbens_basal_ganglia)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region nucleus_accumbens_basal_ganglia`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3966 | 0.1957 | 0.2053 | 0.87 | 0.79 | 0.97 | 8.24e-03 | logit_matched_standard |
| sQTL | all_modules | 3966 | 0.1957 | 0.2053 | 0.87 | 0.78 | 0.97 | 1.18e-02 | logit_matched_constraint |
| eQTL | all_modules | 4318 | 0.465 | 0.4924 | 0.88 | 0.81 | 0.94 | 3.85e-04 | logit_matched_standard |
| eQTL | all_modules | 4318 | 0.465 | 0.4924 | 0.93 | 0.86 | 1.0 | 6.53e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
