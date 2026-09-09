# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6218 | 0.1584 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2531 | 0.1612 | 0.1565 | 1.06 | 0.92 | 1.22 | 4.18e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 566 | 0.1855 | 0.1557 | 1.22 | 0.97 | 1.54 | 8.55e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 1965 | 0.1542 | 0.1604 | 0.98 | 0.84 | 1.14 | 8.29e-01 | logit_matched_standard |
| eQTL | all_modules | 7672 | 0.3693 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3218 | 0.3533 | 0.3808 | 0.9 | 0.82 | 0.98 | 2.29e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 652 | 0.3604 | 0.3701 | 0.97 | 0.82 | 1.15 | 7.30e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 2566 | 0.3515 | 0.3782 | 0.9 | 0.81 | 0.99 | 2.95e-02 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
