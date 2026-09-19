# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13246 | 0.1827 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 755 | 0.1695 | 0.1835 | 0.93 | 0.86 | 1.0 | 4.15e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 755 | 0.1695 | 0.1835 | 0.93 | 0.86 | 1.0 | 4.15e-02 | ols_rankint_matched_standard |
| eQTL | all_modules | 17008 | 0.3873 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 953 | 0.3431 | 0.39 | 0.87 | 0.82 | 0.93 | 3.36e-05 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 953 | 0.3431 | 0.39 | 0.87 | 0.82 | 0.93 | 3.36e-05 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
