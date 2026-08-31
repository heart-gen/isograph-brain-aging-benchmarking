# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6283 | 0.1577 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2898 | 0.1618 | 0.1542 | 0.99 | 0.95 | 1.04 | 8.16e-01 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 713 | 0.1823 | 0.1546 | 1.01 | 0.94 | 1.09 | 7.56e-01 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 2185 | 0.1551 | 0.1591 | 0.99 | 0.94 | 1.04 | 6.52e-01 | ols_rankint_matched_standard |
| eQTL | all_modules | 7736 | 0.3692 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3626 | 0.3593 | 0.3779 | 0.96 | 0.92 | 1.0 | 6.13e-02 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 845 | 0.3787 | 0.368 | 1.01 | 0.94 | 1.09 | 6.88e-01 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 2781 | 0.3535 | 0.378 | 0.95 | 0.91 | 0.99 | 2.73e-02 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
