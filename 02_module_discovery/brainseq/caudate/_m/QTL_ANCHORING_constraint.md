# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-aging/caudate)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region caudate`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4860 | 0.2111 | 0.2067 | 1.03 | 0.94 | 1.14 | 5.02e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 579 | 0.2176 | 0.2081 | 0.95 | 0.77 | 1.18 | 6.60e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 579 | 0.2176 | 0.2081 | 0.95 | 0.77 | 1.18 | 6.60e-01 | logit_matched_standard |
| sQTL | all_modules | 4860 | 0.2111 | 0.2067 | 1.0 | 0.9 | 1.11 | 9.45e-01 | logit_matched_constraint |
| sQTL | pheno_sig_modules | 579 | 0.2176 | 0.2081 | 0.94 | 0.76 | 1.16 | 5.67e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 579 | 0.2176 | 0.2081 | 0.94 | 0.76 | 1.16 | 5.67e-01 | logit_matched_constraint |
| eQTL | all_modules | 5375 | 0.4793 | 0.5225 | 0.83 | 0.77 | 0.89 | 8.02e-08 | logit_matched_standard |
| eQTL | pheno_sig_modules | 630 | 0.4921 | 0.5059 | 0.93 | 0.79 | 1.09 | 3.84e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 630 | 0.4921 | 0.5059 | 0.93 | 0.79 | 1.09 | 3.84e-01 | logit_matched_standard |
| eQTL | all_modules | 5375 | 0.4793 | 0.5225 | 0.87 | 0.8 | 0.94 | 2.54e-04 | logit_matched_constraint |
| eQTL | pheno_sig_modules | 630 | 0.4921 | 0.5059 | 0.94 | 0.8 | 1.11 | 4.91e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 630 | 0.4921 | 0.5059 | 0.94 | 0.8 | 1.11 | 4.91e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
