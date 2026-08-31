# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2341 | 0.217 | 0.2257 | 0.9 | 0.85 | 0.96 | 5.97e-04 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 207 | 0.2512 | 0.2238 | 1.03 | 0.87 | 1.22 | 7.61e-01 | poisson_matched_standard |
| sQTL | go_invisible_modules | 136 | 0.2868 | 0.2236 | 1.07 | 0.88 | 1.3 | 4.83e-01 | poisson_matched_standard |
| sQTL | go_visible_modules | 71 | 0.1831 | 0.2245 | 0.91 | 0.64 | 1.28 | 5.79e-01 | poisson_matched_standard |
| eQTL | all_modules | 3350 | 0.463 | 0.5326 | 0.76 | 0.71 | 0.81 | 2.25e-17 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 285 | 0.407 | 0.5222 | 0.6 | 0.48 | 0.75 | 7.64e-06 | poisson_matched_standard |
| eQTL | go_invisible_modules | 165 | 0.5212 | 0.5205 | 0.78 | 0.6 | 1.01 | 6.00e-02 | poisson_matched_standard |
| eQTL | go_visible_modules | 120 | 0.25 | 0.5222 | 0.34 | 0.22 | 0.54 | 5.67e-06 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
