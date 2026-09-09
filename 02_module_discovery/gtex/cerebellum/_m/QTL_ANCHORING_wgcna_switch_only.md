# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 112 | 0.375 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 112 | 0.375 | nan | nan | nan | nan | NA | insufficient |
| sQTL | go_invisible_modules | 94 | 0.4362 | 0.0556 | 7.22 | 0.78 | 66.77 | 8.17e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 18 | 0.0556 | 0.4362 | 0.14 | 0.01 | 1.28 | 8.17e-02 | logit_matched_standard |
| eQTL | all_modules | 132 | 0.5682 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 132 | 0.5682 | nan | nan | nan | nan | NA | insufficient |
| eQTL | go_invisible_modules | 95 | 0.6947 | 0.2432 | 6.12 | 2.38 | 15.74 | 1.72e-04 | logit_matched_standard |
| eQTL | go_visible_modules | 37 | 0.2432 | 0.6947 | 0.16 | 0.06 | 0.42 | 1.72e-04 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
