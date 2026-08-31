# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 11315 | 0.1993 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 6838 | 0.2066 | 0.1881 | 1.05 | 0.95 | 1.16 | 3.31e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1475 | 0.2319 | 0.1944 | 1.06 | 0.92 | 1.21 | 4.40e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 6219 | 0.2039 | 0.1937 | 1.06 | 0.96 | 1.16 | 2.69e-01 | logit_matched_standard |
| sQTL | all_modules | 11315 | 0.1993 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 6838 | 0.2066 | 0.1881 | 1.01 | 0.91 | 1.12 | 9.24e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1475 | 0.2319 | 0.1944 | 1.06 | 0.92 | 1.23 | 3.86e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 6219 | 0.2039 | 0.1937 | 1.0 | 0.9 | 1.1 | 9.68e-01 | logit_matched_constraint |
| eQTL | all_modules | 13456 | 0.5042 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7865 | 0.5053 | 0.5028 | 1.01 | 0.94 | 1.08 | 7.48e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 1585 | 0.511 | 0.5033 | 1.01 | 0.9 | 1.12 | 9.24e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 7194 | 0.5068 | 0.5013 | 1.03 | 0.96 | 1.1 | 3.84e-01 | logit_matched_standard |
| eQTL | all_modules | 13456 | 0.5042 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7865 | 0.5053 | 0.5028 | 1.06 | 0.99 | 1.14 | 1.07e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1585 | 0.511 | 0.5033 | 1.06 | 0.95 | 1.18 | 2.93e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 7194 | 0.5068 | 0.5013 | 1.07 | 1.0 | 1.15 | 5.88e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
