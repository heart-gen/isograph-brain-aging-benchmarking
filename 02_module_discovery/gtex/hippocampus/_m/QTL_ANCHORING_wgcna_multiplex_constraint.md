# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 10967 | 0.1519 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3654 | 0.1188 | 0.1685 | 0.73 | 0.64 | 0.82 | 2.74e-07 | logit_matched_standard |
| sQTL | go_invisible_modules | 180 | 0.1 | 0.1528 | 0.64 | 0.39 | 1.06 | 8.58e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 3535 | 0.1191 | 0.1675 | 0.73 | 0.65 | 0.83 | 7.48e-07 | logit_matched_standard |
| sQTL | all_modules | 10967 | 0.1519 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3654 | 0.1188 | 0.1685 | 0.74 | 0.66 | 0.84 | 2.83e-06 | logit_matched_constraint |
| sQTL | go_invisible_modules | 180 | 0.1 | 0.1528 | 0.68 | 0.41 | 1.14 | 1.44e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 3535 | 0.1191 | 0.1675 | 0.75 | 0.66 | 0.85 | 4.99e-06 | logit_matched_constraint |
| eQTL | all_modules | 13196 | 0.3634 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4673 | 0.3084 | 0.3936 | 0.69 | 0.64 | 0.75 | 2.49e-21 | logit_matched_standard |
| eQTL | go_invisible_modules | 223 | 0.2556 | 0.3653 | 0.59 | 0.44 | 0.8 | 7.10e-04 | logit_matched_standard |
| eQTL | go_visible_modules | 4515 | 0.3103 | 0.3911 | 0.71 | 0.65 | 0.76 | 5.85e-19 | logit_matched_standard |
| eQTL | all_modules | 13196 | 0.3634 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4673 | 0.3084 | 0.3936 | 0.71 | 0.66 | 0.77 | 3.25e-18 | logit_matched_constraint |
| eQTL | go_invisible_modules | 223 | 0.2556 | 0.3653 | 0.63 | 0.46 | 0.86 | 3.56e-03 | logit_matched_constraint |
| eQTL | go_visible_modules | 4515 | 0.3103 | 0.3911 | 0.72 | 0.67 | 0.78 | 2.20e-16 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
