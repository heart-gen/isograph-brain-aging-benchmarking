# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Threshold-free sensitivity (OLS: rank-inverse-normal of -log10 permutation p ~ module membership + covariates). Uses the same gene-level permutation statistic the sGene/eGene call thresholds, so it adds precision without adding data.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 12905 | 0.1676 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4267 | 0.1362 | 0.1831 | 0.94 | 0.9 | 0.97 | 3.06e-04 | ols_rankint_matched_standard |
| sQTL | go_invisible_modules | 217 | 0.1152 | 0.1685 | 0.86 | 0.75 | 0.98 | 2.26e-02 | ols_rankint_matched_standard |
| sQTL | go_visible_modules | 4117 | 0.1368 | 0.1821 | 0.94 | 0.91 | 0.98 | 1.17e-03 | ols_rankint_matched_standard |
| eQTL | all_modules | 17207 | 0.3904 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 6012 | 0.3395 | 0.4178 | 0.85 | 0.82 | 0.88 | 1.85e-24 | ols_rankint_matched_standard |
| eQTL | go_invisible_modules | 287 | 0.2683 | 0.3925 | 0.78 | 0.7 | 0.88 | 4.24e-05 | ols_rankint_matched_standard |
| eQTL | go_visible_modules | 5796 | 0.3421 | 0.415 | 0.86 | 0.83 | 0.89 | 3.94e-21 | ols_rankint_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
