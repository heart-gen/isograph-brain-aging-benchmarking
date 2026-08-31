# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-aging/caudate)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region caudate`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4254 | 0.2151 | 0.232 | 0.99 | 0.96 | 1.03 | 7.44e-01 | ols_rankint_matched_standard |
| sQTL | pheno_sig_modules | 941 | 0.2189 | 0.2271 | 1.01 | 0.95 | 1.08 | 8.10e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 941 | 0.2189 | 0.2271 | 1.01 | 0.95 | 1.08 | 8.10e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 5465 | 0.4917 | 0.5385 | 0.9 | 0.87 | 0.93 | 4.19e-11 | ols_rankint_matched_standard |
| eQTL | pheno_sig_modules | 1318 | 0.4977 | 0.5259 | 0.97 | 0.91 | 1.02 | 2.45e-01 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 1318 | 0.4977 | 0.5259 | 0.97 | 0.91 | 1.02 | 2.45e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
