# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 11034 | 0.1092 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4723 | 0.1107 | 0.1081 | 1.02 | 0.9 | 1.16 | 7.44e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 733 | 0.1173 | 0.1086 | 0.9 | 0.7 | 1.14 | 3.81e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 4378 | 0.1083 | 0.1098 | 1.01 | 0.89 | 1.15 | 8.35e-01 | logit_matched_standard |
| sQTL | all_modules | 11034 | 0.1092 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4723 | 0.1107 | 0.1081 | 0.96 | 0.85 | 1.09 | 5.35e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 733 | 0.1173 | 0.1086 | 0.84 | 0.66 | 1.08 | 1.73e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 4378 | 0.1083 | 0.1098 | 0.96 | 0.84 | 1.09 | 5.36e-01 | logit_matched_constraint |
| eQTL | all_modules | 13663 | 0.2543 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5830 | 0.2365 | 0.2676 | 0.85 | 0.78 | 0.92 | 4.15e-05 | logit_matched_standard |
| eQTL | go_invisible_modules | 833 | 0.2401 | 0.2553 | 0.9 | 0.76 | 1.06 | 1.98e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 5444 | 0.2333 | 0.2683 | 0.83 | 0.77 | 0.9 | 8.06e-06 | logit_matched_standard |
| eQTL | all_modules | 13663 | 0.2543 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5830 | 0.2365 | 0.2676 | 0.83 | 0.76 | 0.9 | 3.54e-06 | logit_matched_constraint |
| eQTL | go_invisible_modules | 833 | 0.2401 | 0.2553 | 0.87 | 0.73 | 1.02 | 9.00e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 5444 | 0.2333 | 0.2683 | 0.82 | 0.76 | 0.89 | 1.47e-06 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
