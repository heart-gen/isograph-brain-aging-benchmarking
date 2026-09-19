# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 11252 | 0.1668 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 691 | 0.165 | 0.1669 | 0.92 | 0.74 | 1.14 | 4.59e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 691 | 0.165 | 0.1669 | 0.92 | 0.74 | 1.14 | 4.59e-01 | logit_matched_standard |
| sQTL | all_modules | 11252 | 0.1668 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 691 | 0.165 | 0.1669 | 1.0 | 0.8 | 1.25 | 9.95e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 691 | 0.165 | 0.1669 | 1.0 | 0.8 | 1.25 | 9.95e-01 | logit_matched_constraint |
| eQTL | all_modules | 13105 | 0.3599 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 843 | 0.3416 | 0.3612 | 0.88 | 0.76 | 1.01 | 7.76e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 843 | 0.3416 | 0.3612 | 0.88 | 0.76 | 1.01 | 7.76e-02 | logit_matched_standard |
| eQTL | all_modules | 13105 | 0.3599 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 843 | 0.3416 | 0.3612 | 0.88 | 0.76 | 1.02 | 9.08e-02 | logit_matched_constraint |
| eQTL | go_visible_modules | 843 | 0.3416 | 0.3612 | 0.88 | 0.76 | 1.02 | 9.08e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
