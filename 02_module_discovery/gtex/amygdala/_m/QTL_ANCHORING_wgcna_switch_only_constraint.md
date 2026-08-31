# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4238 | 0.1024 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3297 | 0.0986 | 0.1158 | 0.99 | 0.78 | 1.26 | 9.57e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1805 | 0.1014 | 0.1032 | 1.03 | 0.84 | 1.27 | 7.92e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1492 | 0.0952 | 0.1063 | 0.96 | 0.78 | 1.2 | 7.47e-01 | logit_matched_standard |
| sQTL | all_modules | 4238 | 0.1024 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3297 | 0.0986 | 0.1158 | 0.98 | 0.77 | 1.26 | 8.77e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1805 | 0.1014 | 0.1032 | 1.03 | 0.83 | 1.28 | 7.72e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 1492 | 0.0952 | 0.1063 | 0.95 | 0.76 | 1.19 | 6.61e-01 | logit_matched_constraint |
| eQTL | all_modules | 5052 | 0.2381 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3876 | 0.2371 | 0.2415 | 1.0 | 0.85 | 1.16 | 9.50e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 1991 | 0.2396 | 0.2372 | 1.02 | 0.89 | 1.16 | 8.20e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1885 | 0.2345 | 0.2403 | 0.98 | 0.86 | 1.12 | 7.76e-01 | logit_matched_standard |
| eQTL | all_modules | 5052 | 0.2381 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3876 | 0.2371 | 0.2415 | 1.02 | 0.87 | 1.19 | 8.01e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1991 | 0.2396 | 0.2372 | 1.04 | 0.91 | 1.2 | 5.40e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 1885 | 0.2345 | 0.2403 | 0.97 | 0.85 | 1.12 | 6.96e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
