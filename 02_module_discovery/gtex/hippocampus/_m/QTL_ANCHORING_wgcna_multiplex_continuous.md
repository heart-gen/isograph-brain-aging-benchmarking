# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 12773 | 0.1666 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2792 | 0.1357 | 0.1752 | 0.96 | 0.92 | 1.0 | 3.56e-02 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 267 | 0.1723 | 0.1665 | 0.94 | 0.84 | 1.06 | 3.29e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 2601 | 0.1319 | 0.1755 | 0.96 | 0.92 | 1.0 | 4.26e-02 | ols_rankint_matched_standard |
| eQTL | all_modules | 16698 | 0.3894 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4002 | 0.3143 | 0.413 | 0.82 | 0.79 | 0.85 | 3.62e-27 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 339 | 0.3097 | 0.391 | 0.86 | 0.77 | 0.96 | 6.11e-03 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 3742 | 0.3135 | 0.4113 | 0.83 | 0.8 | 0.86 | 3.14e-25 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
