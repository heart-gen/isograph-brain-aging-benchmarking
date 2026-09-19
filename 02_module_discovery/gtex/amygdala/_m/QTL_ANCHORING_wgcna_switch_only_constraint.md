# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4726 | 0.1185 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2471 | 0.1141 | 0.1233 | 0.91 | 0.75 | 1.09 | 2.87e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1406 | 0.1117 | 0.1214 | 0.92 | 0.76 | 1.13 | 4.48e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1065 | 0.1174 | 0.1188 | 0.95 | 0.77 | 1.18 | 6.57e-01 | logit_matched_standard |
| sQTL | all_modules | 4726 | 0.1185 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2471 | 0.1141 | 0.1233 | 0.91 | 0.76 | 1.1 | 3.20e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1406 | 0.1117 | 0.1214 | 0.94 | 0.77 | 1.16 | 5.72e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1065 | 0.1174 | 0.1188 | 0.94 | 0.75 | 1.17 | 5.70e-01 | logit_matched_constraint |
| eQTL | all_modules | 5092 | 0.2545 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 2652 | 0.2579 | 0.2508 | 1.03 | 0.91 | 1.17 | 6.72e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 1518 | 0.2576 | 0.2532 | 1.02 | 0.89 | 1.17 | 7.87e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1134 | 0.2584 | 0.2534 | 1.02 | 0.87 | 1.18 | 8.33e-01 | logit_matched_standard |
| eQTL | all_modules | 5092 | 0.2545 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 2652 | 0.2579 | 0.2508 | 1.0 | 0.88 | 1.14 | 9.42e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1518 | 0.2576 | 0.2532 | 1.0 | 0.87 | 1.15 | 9.73e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1134 | 0.2584 | 0.2534 | 1.0 | 0.86 | 1.17 | 9.60e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
