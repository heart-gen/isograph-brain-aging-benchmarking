# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-aging/caudate)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region caudate`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3714 | 0.2046 | 0.2136 | 0.97 | 0.87 | 1.07 | 5.22e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 809 | 0.2064 | 0.211 | 0.98 | 0.82 | 1.18 | 8.33e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 809 | 0.2064 | 0.211 | 0.98 | 0.82 | 1.18 | 8.33e-01 | logit_matched_standard |
| sQTL | all_modules | 3714 | 0.2046 | 0.2136 | 0.97 | 0.88 | 1.08 | 6.14e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 809 | 0.2064 | 0.211 | 0.97 | 0.8 | 1.17 | 7.44e-01 | logit_matched_constraint |
| sQTL | go_invisible_modules | 809 | 0.2064 | 0.211 | 0.97 | 0.8 | 1.17 | 7.44e-01 | logit_matched_constraint |
| eQTL | all_modules | 4203 | 0.4916 | 0.5117 | 0.92 | 0.85 | 0.99 | 1.93e-02 | logit_matched_standard |
| eQTL | pheno_sig_modules | 887 | 0.5152 | 0.5045 | 1.04 | 0.9 | 1.19 | 5.95e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 887 | 0.5152 | 0.5045 | 1.04 | 0.9 | 1.19 | 5.95e-01 | logit_matched_standard |
| eQTL | all_modules | 4203 | 0.4916 | 0.5117 | 0.97 | 0.9 | 1.05 | 4.34e-01 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 887 | 0.5152 | 0.5045 | 1.06 | 0.92 | 1.22 | 3.92e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 887 | 0.5152 | 0.5045 | 1.06 | 0.92 | 1.22 | 3.92e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
