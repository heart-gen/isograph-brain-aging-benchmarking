# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 734 | 0.3283 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 532 | 0.3477 | 0.2772 | 0.88 | 0.59 | 1.32 | 5.45e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 236 | 0.3686 | 0.3092 | 0.95 | 0.66 | 1.36 | 7.60e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 296 | 0.3311 | 0.3265 | 0.96 | 0.68 | 1.36 | 8.22e-01 | logit_matched_standard |
| sQTL | all_modules | 734 | 0.3283 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 532 | 0.3477 | 0.2772 | 1.13 | 0.74 | 1.73 | 5.80e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 236 | 0.3686 | 0.3092 | 1.07 | 0.73 | 1.56 | 7.41e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 296 | 0.3311 | 0.3265 | 1.03 | 0.72 | 1.47 | 8.74e-01 | logit_matched_constraint |
| eQTL | all_modules | 855 | 0.5801 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 598 | 0.5936 | 0.5486 | 1.14 | 0.85 | 1.55 | 3.80e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 274 | 0.5876 | 0.5766 | 1.02 | 0.76 | 1.37 | 8.91e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 324 | 0.5988 | 0.5687 | 1.1 | 0.83 | 1.47 | 4.92e-01 | logit_matched_standard |
| eQTL | all_modules | 855 | 0.5801 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 598 | 0.5936 | 0.5486 | 1.26 | 0.92 | 1.73 | 1.55e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 274 | 0.5876 | 0.5766 | 1.11 | 0.82 | 1.51 | 4.95e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 324 | 0.5988 | 0.5687 | 1.1 | 0.82 | 1.49 | 5.11e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
