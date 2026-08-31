# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Nucleus_accumbens_basal_ganglia xQTL (gtex-aging/nucleus_accumbens_basal_ganglia)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region nucleus_accumbens_basal_ganglia`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4731 | 0.2228 | 0.2187 | 0.96 | 0.93 | 1.0 | 4.19e-02 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 104 | 0.0577 | 0.2213 | 0.77 | 0.64 | 0.93 | 6.51e-03 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 104 | 0.0577 | 0.2213 | 0.77 | 0.64 | 0.93 | 6.51e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 5690 | 0.49 | 0.52 | 0.9 | 0.87 | 0.93 | 2.42e-11 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 194 | 0.3557 | 0.5124 | 0.66 | 0.57 | 0.76 | 8.79e-09 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 194 | 0.3557 | 0.5124 | 0.66 | 0.57 | 0.76 | 8.79e-09 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
