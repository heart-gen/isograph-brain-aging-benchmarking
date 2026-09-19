# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5649 | 0.2261 | 0.2236 | 0.94 | 0.9 | 0.99 | 1.05e-02 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 1618 | 0.2831 | 0.2167 | 1.19 | 1.12 | 1.27 | 1.93e-08 | poisson_matched_standard |
| sQTL | go_invisible_modules | 1297 | 0.2791 | 0.2189 | 1.1 | 1.03 | 1.18 | 4.80e-03 | poisson_matched_standard |
| sQTL | go_visible_modules | 321 | 0.2991 | 0.2228 | 1.55 | 1.37 | 1.76 | 4.09e-12 | poisson_matched_standard |
| eQTL | all_modules | 6983 | 0.4875 | 0.5476 | 0.79 | 0.75 | 0.82 | 2.58e-23 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 1809 | 0.5412 | 0.5232 | 1.06 | 0.99 | 1.14 | 1.09e-01 | poisson_matched_standard |
| eQTL | go_invisible_modules | 1427 | 0.5501 | 0.5229 | 1.1 | 1.01 | 1.19 | 2.35e-02 | poisson_matched_standard |
| eQTL | go_visible_modules | 382 | 0.5079 | 0.5253 | 0.93 | 0.79 | 1.09 | 3.57e-01 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
