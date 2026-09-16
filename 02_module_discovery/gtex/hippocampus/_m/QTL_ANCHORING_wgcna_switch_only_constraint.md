# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5108 | 0.1558 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2186 | 0.1889 | 0.1311 | 1.32 | 1.13 | 1.54 | 6.07e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 1950 | 0.1897 | 0.1349 | 1.3 | 1.11 | 1.52 | 1.36e-03 | logit_matched_standard |
| sQTL | go_visible_modules | 236 | 0.1822 | 0.1546 | 1.12 | 0.79 | 1.59 | 5.31e-01 | logit_matched_standard |
| sQTL | all_modules | 5108 | 0.1558 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2186 | 0.1889 | 0.1311 | 1.25 | 1.07 | 1.47 | 5.59e-03 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1950 | 0.1897 | 0.1349 | 1.24 | 1.06 | 1.46 | 8.37e-03 | logit_matched_constraint |
| sQTL | go_visible_modules | 236 | 0.1822 | 0.1546 | 1.07 | 0.75 | 1.53 | 7.01e-01 | logit_matched_constraint |
| eQTL | all_modules | 5493 | 0.3592 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 2318 | 0.3745 | 0.348 | 1.11 | 0.99 | 1.24 | 7.17e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 2072 | 0.3764 | 0.3487 | 1.11 | 0.99 | 1.25 | 7.11e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 246 | 0.3577 | 0.3593 | 1.01 | 0.77 | 1.32 | 9.36e-01 | logit_matched_standard |
| eQTL | all_modules | 5493 | 0.3592 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 2318 | 0.3745 | 0.348 | 1.05 | 0.93 | 1.17 | 4.32e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 2072 | 0.3764 | 0.3487 | 1.05 | 0.94 | 1.18 | 3.96e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 246 | 0.3577 | 0.3593 | 0.99 | 0.75 | 1.29 | 9.14e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
