# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cerebellar_Hemisphere xQTL (gtex-aging/cerebellar_hemisphere)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellar_hemisphere --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 10302 | 0.2831 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 114 | 0.1491 | 0.2846 | 0.74 | 0.43 | 1.28 | 2.80e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 114 | 0.1491 | 0.2846 | 0.74 | 0.43 | 1.28 | 2.80e-01 | logit_matched_standard |
| sQTL | all_modules | 10302 | 0.2831 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 114 | 0.1491 | 0.2846 | 0.74 | 0.43 | 1.3 | 2.99e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 114 | 0.1491 | 0.2846 | 0.74 | 0.43 | 1.3 | 2.99e-01 | logit_matched_constraint |
| eQTL | all_modules | 12073 | 0.5898 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 172 | 0.4186 | 0.5923 | 0.48 | 0.35 | 0.65 | 2.54e-06 | logit_matched_standard |
| eQTL | go_visible_modules | 172 | 0.4186 | 0.5923 | 0.48 | 0.35 | 0.65 | 2.54e-06 | logit_matched_standard |
| eQTL | all_modules | 12073 | 0.5898 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 172 | 0.4186 | 0.5923 | 0.53 | 0.39 | 0.72 | 6.90e-05 | logit_matched_constraint |
| eQTL | go_visible_modules | 172 | 0.4186 | 0.5923 | 0.53 | 0.39 | 0.72 | 6.90e-05 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
