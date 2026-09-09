# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 9966 | 0.2673 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3585 | 0.3026 | 0.2475 | 1.19 | 1.08 | 1.31 | 7.05e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 299 | 0.3645 | 0.2643 | 1.46 | 1.13 | 1.9 | 4.30e-03 | logit_matched_standard |
| sQTL | go_visible_modules | 3296 | 0.297 | 0.2526 | 1.13 | 1.02 | 1.25 | 1.88e-02 | logit_matched_standard |
| eQTL | all_modules | 12763 | 0.5771 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4540 | 0.5921 | 0.5689 | 1.1 | 1.02 | 1.19 | 1.06e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 385 | 0.7273 | 0.5725 | 1.97 | 1.57 | 2.48 | 5.10e-09 | logit_matched_standard |
| eQTL | go_visible_modules | 4165 | 0.5798 | 0.5758 | 1.02 | 0.94 | 1.1 | 6.37e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
