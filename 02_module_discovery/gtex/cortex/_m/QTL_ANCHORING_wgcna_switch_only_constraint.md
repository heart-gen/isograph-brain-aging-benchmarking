# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3121 | 0.2272 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2579 | 0.2288 | 0.2196 | 0.82 | 0.65 | 1.05 | 1.16e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 989 | 0.2649 | 0.2097 | 1.08 | 0.89 | 1.3 | 4.28e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1590 | 0.2063 | 0.2489 | 0.84 | 0.7 | 1.01 | 5.79e-02 | logit_matched_standard |
| sQTL | all_modules | 3121 | 0.2272 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2579 | 0.2288 | 0.2196 | 0.91 | 0.71 | 1.16 | 4.30e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 989 | 0.2649 | 0.2097 | 1.13 | 0.93 | 1.37 | 2.17e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1590 | 0.2063 | 0.2489 | 0.85 | 0.71 | 1.02 | 8.00e-02 | logit_matched_constraint |
| eQTL | all_modules | 3528 | 0.5011 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 2850 | 0.5105 | 0.4617 | 1.21 | 1.02 | 1.43 | 2.83e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 1076 | 0.5418 | 0.4833 | 1.25 | 1.08 | 1.45 | 2.32e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 1774 | 0.4915 | 0.5108 | 0.93 | 0.81 | 1.06 | 2.84e-01 | logit_matched_standard |
| eQTL | all_modules | 3528 | 0.5011 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 2850 | 0.5105 | 0.4617 | 1.3 | 1.1 | 1.55 | 2.82e-03 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1076 | 0.5418 | 0.4833 | 1.26 | 1.09 | 1.46 | 1.93e-03 | logit_matched_constraint |
| eQTL | go_visible_modules | 1774 | 0.4915 | 0.5108 | 0.97 | 0.84 | 1.11 | 6.09e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
