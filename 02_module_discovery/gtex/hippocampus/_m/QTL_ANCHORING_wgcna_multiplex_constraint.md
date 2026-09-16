# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 10918 | 0.1516 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2381 | 0.1163 | 0.1614 | 0.72 | 0.63 | 0.84 | 9.87e-06 | logit_matched_standard |
| sQTL | go_invisible_modules | 225 | 0.1689 | 0.1512 | 0.93 | 0.64 | 1.34 | 6.86e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 2226 | 0.1114 | 0.1619 | 0.69 | 0.6 | 0.8 | 1.39e-06 | logit_matched_standard |
| sQTL | all_modules | 10918 | 0.1516 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2381 | 0.1163 | 0.1614 | 0.75 | 0.64 | 0.86 | 7.11e-05 | logit_matched_constraint |
| sQTL | go_invisible_modules | 225 | 0.1689 | 0.1512 | 0.93 | 0.64 | 1.36 | 7.23e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 2226 | 0.1114 | 0.1619 | 0.71 | 0.61 | 0.83 | 1.24e-05 | logit_matched_constraint |
| eQTL | all_modules | 13051 | 0.3636 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3140 | 0.2806 | 0.3899 | 0.61 | 0.56 | 0.67 | 1.57e-27 | logit_matched_standard |
| eQTL | go_invisible_modules | 264 | 0.2803 | 0.3653 | 0.67 | 0.51 | 0.88 | 3.83e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 2948 | 0.2795 | 0.3881 | 0.61 | 0.56 | 0.67 | 4.33e-26 | logit_matched_standard |
| eQTL | all_modules | 13051 | 0.3636 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3140 | 0.2806 | 0.3899 | 0.62 | 0.56 | 0.67 | 3.65e-26 | logit_matched_constraint |
| eQTL | go_invisible_modules | 264 | 0.2803 | 0.3653 | 0.7 | 0.53 | 0.92 | 9.73e-03 | logit_matched_constraint |
| eQTL | go_visible_modules | 2948 | 0.2795 | 0.3881 | 0.62 | 0.56 | 0.68 | 8.93e-25 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
