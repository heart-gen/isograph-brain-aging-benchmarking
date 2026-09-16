# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 11298 | 0.199 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7087 | 0.2119 | 0.1772 | 1.09 | 0.99 | 1.21 | 9.09e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 2748 | 0.1921 | 0.2012 | 0.84 | 0.75 | 0.94 | 2.26e-03 | logit_matched_standard |
| sQTL | go_visible_modules | 5668 | 0.2239 | 0.1739 | 1.25 | 1.13 | 1.38 | 7.20e-06 | logit_matched_standard |
| sQTL | all_modules | 11298 | 0.199 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7087 | 0.2119 | 0.1772 | 1.04 | 0.94 | 1.16 | 4.42e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 2748 | 0.1921 | 0.2012 | 0.85 | 0.75 | 0.95 | 5.52e-03 | logit_matched_constraint |
| sQTL | go_visible_modules | 5668 | 0.2239 | 0.1739 | 1.18 | 1.06 | 1.3 | 1.76e-03 | logit_matched_constraint |
| eQTL | all_modules | 13367 | 0.5053 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7998 | 0.5134 | 0.4934 | 1.08 | 1.01 | 1.16 | 3.33e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 2968 | 0.4912 | 0.5094 | 0.92 | 0.84 | 1.0 | 3.95e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 6430 | 0.5202 | 0.4916 | 1.12 | 1.05 | 1.2 | 1.19e-03 | logit_matched_standard |
| eQTL | all_modules | 13367 | 0.5053 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7998 | 0.5134 | 0.4934 | 1.15 | 1.07 | 1.24 | 1.78e-04 | logit_matched_constraint |
| eQTL | go_invisible_modules | 2968 | 0.4912 | 0.5094 | 0.98 | 0.9 | 1.07 | 7.18e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 6430 | 0.5202 | 0.4916 | 1.16 | 1.08 | 1.25 | 4.99e-05 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
