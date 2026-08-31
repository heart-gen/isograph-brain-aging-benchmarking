# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 8830 | 0.2613 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 6823 | 0.2555 | 0.281 | 0.86 | 0.76 | 0.97 | 1.49e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 227 | 0.3612 | 0.2586 | 1.36 | 1.01 | 1.83 | 4.43e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 6634 | 0.2532 | 0.2855 | 0.84 | 0.75 | 0.95 | 4.67e-03 | logit_matched_standard |
| sQTL | all_modules | 8830 | 0.2613 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 6823 | 0.2555 | 0.281 | 0.85 | 0.75 | 0.96 | 1.18e-02 | logit_matched_constraint |
| sQTL | go_invisible_modules | 227 | 0.3612 | 0.2586 | 1.16 | 0.85 | 1.57 | 3.50e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 6634 | 0.2532 | 0.2855 | 0.85 | 0.76 | 0.96 | 1.05e-02 | logit_matched_constraint |
| eQTL | all_modules | 10436 | 0.5728 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7752 | 0.5787 | 0.5559 | 1.08 | 0.99 | 1.18 | 8.17e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 247 | 0.7409 | 0.5688 | 2.12 | 1.59 | 2.83 | 3.17e-07 | logit_matched_standard |
| eQTL | go_visible_modules | 7544 | 0.5741 | 0.5695 | 1.01 | 0.92 | 1.1 | 8.79e-01 | logit_matched_standard |
| eQTL | all_modules | 10436 | 0.5728 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7752 | 0.5787 | 0.5559 | 0.98 | 0.89 | 1.08 | 6.64e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 247 | 0.7409 | 0.5688 | 1.81 | 1.35 | 2.43 | 6.89e-05 | logit_matched_constraint |
| eQTL | go_visible_modules | 7544 | 0.5741 | 0.5695 | 0.93 | 0.85 | 1.01 | 9.96e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
