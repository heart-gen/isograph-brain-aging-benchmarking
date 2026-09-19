# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (brainseq-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3976 | 0.1461 | 0.1775 | 0.76 | 0.68 | 0.85 | 6.34e-07 | logit_matched_standard |
| sQTL | pheno_sig_modules | 74 | 0.0946 | 0.1684 | 0.45 | 0.2 | 0.99 | 4.65e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 74 | 0.0946 | 0.1684 | 0.45 | 0.2 | 0.99 | 4.65e-02 | logit_matched_standard |
| eQTL | all_modules | 4588 | 0.3354 | 0.3975 | 0.74 | 0.69 | 0.8 | 4.16e-16 | logit_matched_standard |
| eQTL | pheno_sig_modules | 76 | 0.3421 | 0.3818 | 0.78 | 0.49 | 1.26 | 3.19e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 76 | 0.3421 | 0.3818 | 0.78 | 0.49 | 1.26 | 3.19e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
