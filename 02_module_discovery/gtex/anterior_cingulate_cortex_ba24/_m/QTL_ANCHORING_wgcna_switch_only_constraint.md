# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5575 | 0.143 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3749 | 0.1483 | 0.132 | 1.0 | 0.85 | 1.19 | 9.55e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 79 | 0.2025 | 0.1421 | 1.44 | 0.81 | 2.58 | 2.18e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 3670 | 0.1471 | 0.1349 | 0.98 | 0.83 | 1.15 | 7.90e-01 | logit_matched_standard |
| sQTL | all_modules | 5575 | 0.143 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3749 | 0.1483 | 0.132 | 0.99 | 0.84 | 1.18 | 9.49e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 79 | 0.2025 | 0.1421 | 1.3 | 0.72 | 2.35 | 3.85e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 3670 | 0.1471 | 0.1349 | 0.97 | 0.82 | 1.16 | 7.65e-01 | logit_matched_constraint |
| eQTL | all_modules | 6341 | 0.3574 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4241 | 0.3577 | 0.3567 | 1.01 | 0.9 | 1.12 | 9.15e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 91 | 0.3846 | 0.357 | 1.09 | 0.71 | 1.68 | 6.84e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 4150 | 0.3571 | 0.3578 | 1.0 | 0.9 | 1.11 | 9.96e-01 | logit_matched_standard |
| eQTL | all_modules | 6341 | 0.3574 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4241 | 0.3577 | 0.3567 | 0.99 | 0.88 | 1.11 | 8.49e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 91 | 0.3846 | 0.357 | 1.08 | 0.7 | 1.68 | 7.28e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 4150 | 0.3571 | 0.3578 | 0.98 | 0.88 | 1.1 | 7.83e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
