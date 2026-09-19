# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (brainseq-aging/hippocampus)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3976 | 0.1461 | 0.1775 | 0.82 | 0.77 | 0.87 | 7.05e-11 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 74 | 0.0946 | 0.1684 | 0.71 | 0.47 | 1.06 | 9.13e-02 | poisson_matched_standard |
| sQTL | go_invisible_modules | 74 | 0.0946 | 0.1684 | 0.71 | 0.47 | 1.06 | 9.13e-02 | poisson_matched_standard |
| eQTL | all_modules | 4588 | 0.3354 | 0.3975 | 0.75 | 0.7 | 0.8 | 6.23e-19 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 76 | 0.3421 | 0.3818 | 0.69 | 0.44 | 1.1 | 1.22e-01 | poisson_matched_standard |
| eQTL | go_invisible_modules | 76 | 0.3421 | 0.3818 | 0.69 | 0.44 | 1.1 | 1.22e-01 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
