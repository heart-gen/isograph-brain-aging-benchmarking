# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13048 | 0.1211 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 8326 | 0.1128 | 0.1357 | 0.76 | 0.68 | 0.85 | 1.20e-06 | logit_matched_standard |
| sQTL | go_invisible_modules | 824 | 0.1262 | 0.1207 | 0.82 | 0.66 | 1.03 | 8.54e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 8131 | 0.112 | 0.1361 | 0.76 | 0.68 | 0.85 | 1.23e-06 | logit_matched_standard |
| eQTL | all_modules | 17753 | 0.283 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 10929 | 0.265 | 0.3118 | 0.79 | 0.74 | 0.84 | 2.93e-12 | logit_matched_standard |
| eQTL | go_invisible_modules | 1034 | 0.2524 | 0.2849 | 0.81 | 0.7 | 0.93 | 3.91e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 10620 | 0.2655 | 0.309 | 0.8 | 0.75 | 0.86 | 1.08e-10 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
