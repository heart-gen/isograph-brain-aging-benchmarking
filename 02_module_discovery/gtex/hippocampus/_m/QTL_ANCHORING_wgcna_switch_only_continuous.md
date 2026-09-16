# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5677 | 0.1649 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2485 | 0.1936 | 0.1425 | 1.08 | 1.02 | 1.13 | 5.13e-03 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 2209 | 0.1951 | 0.1456 | 1.06 | 1.01 | 1.12 | 2.16e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 276 | 0.1812 | 0.164 | 1.08 | 0.96 | 1.22 | 2.12e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 6286 | 0.3721 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 2737 | 0.3906 | 0.3578 | 1.05 | 1.0 | 1.11 | 4.16e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 2435 | 0.3914 | 0.3599 | 1.05 | 0.99 | 1.1 | 8.06e-02 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 302 | 0.3841 | 0.3715 | 1.05 | 0.93 | 1.17 | 4.52e-01 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
