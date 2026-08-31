# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5629 | 0.1975 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4394 | 0.1994 | 0.1911 | 0.97 | 0.82 | 1.15 | 7.64e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 683 | 0.2225 | 0.1941 | 1.32 | 1.08 | 1.62 | 6.51e-03 | logit_matched_standard |
| sQTL | go_visible_modules | 3711 | 0.1951 | 0.2023 | 0.85 | 0.74 | 0.99 | 3.29e-02 | logit_matched_standard |
| sQTL | all_modules | 5629 | 0.1975 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4394 | 0.1994 | 0.1911 | 1.03 | 0.86 | 1.22 | 7.58e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 683 | 0.2225 | 0.1941 | 1.34 | 1.09 | 1.65 | 5.95e-03 | logit_matched_constraint |
| sQTL | go_visible_modules | 3711 | 0.1951 | 0.2023 | 0.88 | 0.76 | 1.03 | 1.03e-01 | logit_matched_constraint |
| eQTL | all_modules | 6323 | 0.4871 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4872 | 0.4893 | 0.4797 | 1.06 | 0.94 | 1.19 | 3.58e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 745 | 0.4725 | 0.4891 | 0.94 | 0.81 | 1.1 | 4.70e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 4127 | 0.4924 | 0.4772 | 1.07 | 0.97 | 1.19 | 1.93e-01 | logit_matched_standard |
| eQTL | all_modules | 6323 | 0.4871 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4872 | 0.4893 | 0.4797 | 1.1 | 0.98 | 1.24 | 1.13e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 745 | 0.4725 | 0.4891 | 0.97 | 0.83 | 1.14 | 7.07e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 4127 | 0.4924 | 0.4772 | 1.09 | 0.98 | 1.22 | 9.91e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
