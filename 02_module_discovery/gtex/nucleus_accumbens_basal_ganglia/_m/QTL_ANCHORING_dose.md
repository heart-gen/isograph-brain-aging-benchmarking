# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Nucleus_accumbens_basal_ganglia xQTL (gtex-aging/nucleus_accumbens_basal_ganglia)

Dose sensitivity (Poisson: number of independent SuSiE credible sets ~ module membership + covariates). Distinguishes a gene with one weak cis signal from one with several independent ones; LD-resolved.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region nucleus_accumbens_basal_ganglia`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4731 | 0.2228 | 0.2187 | 0.89 | 0.85 | 0.94 | 4.25e-06 | poisson_matched_standard |
| sQTL | pheno_sig_modules | 104 | 0.0577 | 0.2213 | 0.25 | 0.13 | 0.46 | 9.76e-06 | poisson_matched_standard |
| sQTL | go_visible_modules | 104 | 0.0577 | 0.2213 | 0.25 | 0.13 | 0.46 | 9.76e-06 | poisson_matched_standard |
| eQTL | all_modules | 5690 | 0.49 | 0.52 | 0.88 | 0.84 | 0.92 | 2.74e-07 | poisson_matched_standard |
| eQTL | pheno_sig_modules | 194 | 0.3557 | 0.5124 | 0.45 | 0.33 | 0.61 | 4.16e-07 | poisson_matched_standard |
| eQTL | go_visible_modules | 194 | 0.3557 | 0.5124 | 0.45 | 0.33 | 0.61 | 4.16e-07 | poisson_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
