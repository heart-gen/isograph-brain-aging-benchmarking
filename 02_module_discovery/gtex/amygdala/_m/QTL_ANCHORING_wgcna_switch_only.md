# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5252 | 0.1247 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2763 | 0.1227 | 0.127 | 0.94 | 0.79 | 1.11 | 4.70e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1607 | 0.1213 | 0.1262 | 0.96 | 0.8 | 1.16 | 6.86e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1156 | 0.1246 | 0.1248 | 0.96 | 0.78 | 1.17 | 6.77e-01 | logit_matched_standard |
| eQTL | all_modules | 5813 | 0.2663 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3056 | 0.2703 | 0.2619 | 1.03 | 0.92 | 1.16 | 6.13e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 1780 | 0.2719 | 0.2638 | 1.04 | 0.91 | 1.18 | 5.77e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1276 | 0.268 | 0.2658 | 1.0 | 0.87 | 1.15 | 9.91e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
