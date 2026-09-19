# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (brainseq-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3651 | 0.1408 | 0.1552 | 0.92 | 0.81 | 1.03 | 1.43e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 70 | 0.0857 | 0.1509 | 0.51 | 0.21 | 1.19 | 1.18e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 70 | 0.0857 | 0.1509 | 0.51 | 0.21 | 1.19 | 1.18e-01 | logit_matched_standard |
| sQTL | all_modules | 3651 | 0.1408 | 0.1552 | 0.94 | 0.83 | 1.07 | 3.51e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 70 | 0.0857 | 0.1509 | 0.52 | 0.22 | 1.24 | 1.40e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 70 | 0.0857 | 0.1509 | 0.52 | 0.22 | 1.24 | 1.40e-01 | logit_matched_constraint |
| eQTL | all_modules | 4136 | 0.3274 | 0.3752 | 0.79 | 0.73 | 0.86 | 7.45e-09 | logit_matched_standard |
| eQTL | pheno_sig_modules | 72 | 0.3333 | 0.3606 | 0.85 | 0.52 | 1.38 | 5.05e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 72 | 0.3333 | 0.3606 | 0.85 | 0.52 | 1.38 | 5.05e-01 | logit_matched_standard |
| eQTL | all_modules | 4136 | 0.3274 | 0.3752 | 0.87 | 0.8 | 0.94 | 9.16e-04 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 72 | 0.3333 | 0.3606 | 0.91 | 0.55 | 1.5 | 7.19e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 72 | 0.3333 | 0.3606 | 0.91 | 0.55 | 1.5 | 7.19e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
