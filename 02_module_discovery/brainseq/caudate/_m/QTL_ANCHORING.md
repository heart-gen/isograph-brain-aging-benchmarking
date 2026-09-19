# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-aging/caudate)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region caudate`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5425 | 0.2199 | 0.2278 | 0.92 | 0.84 | 1.0 | 6.28e-02 | logit_matched_standard |
| sQTL | pheno_sig_modules | 665 | 0.2286 | 0.2244 | 0.9 | 0.74 | 1.09 | 2.60e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 665 | 0.2286 | 0.2244 | 0.9 | 0.74 | 1.09 | 2.60e-01 | logit_matched_standard |
| eQTL | all_modules | 6578 | 0.4837 | 0.5471 | 0.77 | 0.72 | 0.82 | 1.79e-17 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1141 | 0.4566 | 0.5291 | 0.8 | 0.71 | 0.9 | 2.75e-04 | logit_matched_standard |
| eQTL | go_visible_modules | 1141 | 0.4566 | 0.5291 | 0.8 | 0.71 | 0.9 | 2.75e-04 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
