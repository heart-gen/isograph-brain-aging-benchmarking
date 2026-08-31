# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 11347 | 0.1671 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 448 | 0.1585 | 0.1674 | 0.96 | 0.73 | 1.26 | 7.83e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 137 | 0.2336 | 0.1663 | 1.19 | 0.78 | 1.81 | 4.28e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 311 | 0.1254 | 0.1683 | 0.84 | 0.59 | 1.2 | 3.41e-01 | logit_matched_standard |
| sQTL | all_modules | 11347 | 0.1671 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 448 | 0.1585 | 0.1674 | 1.08 | 0.81 | 1.42 | 6.08e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 137 | 0.2336 | 0.1663 | 1.21 | 0.79 | 1.85 | 3.88e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 311 | 0.1254 | 0.1683 | 0.99 | 0.69 | 1.43 | 9.62e-01 | logit_matched_constraint |
| eQTL | all_modules | 13358 | 0.359 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 604 | 0.3162 | 0.361 | 0.76 | 0.64 | 0.91 | 2.38e-03 | logit_matched_standard |
| eQTL | go_invisible_modules | 141 | 0.3262 | 0.3593 | 0.76 | 0.53 | 1.08 | 1.23e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 463 | 0.3132 | 0.3606 | 0.77 | 0.63 | 0.94 | 9.63e-03 | logit_matched_standard |
| eQTL | all_modules | 13358 | 0.359 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 604 | 0.3162 | 0.361 | 0.72 | 0.6 | 0.86 | 3.94e-04 | logit_matched_constraint |
| eQTL | go_invisible_modules | 141 | 0.3262 | 0.3593 | 0.82 | 0.57 | 1.19 | 2.97e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 463 | 0.3132 | 0.3606 | 0.69 | 0.56 | 0.85 | 5.39e-04 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
