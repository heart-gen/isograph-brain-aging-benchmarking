# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4806 | 0.0976 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4272 | 0.0974 | 0.0993 | 0.87 | 0.64 | 1.19 | 3.88e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 142 | 0.0704 | 0.0984 | 0.79 | 0.41 | 1.53 | 4.83e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 4130 | 0.0983 | 0.0932 | 0.94 | 0.7 | 1.25 | 6.60e-01 | logit_matched_standard |
| sQTL | all_modules | 4806 | 0.0976 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4272 | 0.0974 | 0.0993 | 0.91 | 0.66 | 1.25 | 5.61e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 142 | 0.0704 | 0.0984 | 0.71 | 0.36 | 1.38 | 3.14e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 4130 | 0.0983 | 0.0932 | 1.0 | 0.74 | 1.34 | 9.84e-01 | logit_matched_constraint |
| eQTL | all_modules | 5594 | 0.2376 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4905 | 0.2318 | 0.2787 | 0.79 | 0.66 | 0.95 | 1.02e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 149 | 0.2081 | 0.2384 | 0.86 | 0.57 | 1.28 | 4.54e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 4756 | 0.2325 | 0.2661 | 0.84 | 0.71 | 0.99 | 4.25e-02 | logit_matched_standard |
| eQTL | all_modules | 5594 | 0.2376 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4905 | 0.2318 | 0.2787 | 0.86 | 0.71 | 1.03 | 1.01e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 149 | 0.2081 | 0.2384 | 0.87 | 0.58 | 1.31 | 4.96e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 4756 | 0.2325 | 0.2661 | 0.9 | 0.76 | 1.07 | 2.28e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
