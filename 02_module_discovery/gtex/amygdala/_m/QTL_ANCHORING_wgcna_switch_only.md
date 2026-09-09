# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4977 | 0.1127 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3916 | 0.107 | 0.1338 | 0.86 | 0.7 | 1.07 | 1.81e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 2064 | 0.1037 | 0.1191 | 0.84 | 0.7 | 1.02 | 7.62e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 1852 | 0.1107 | 0.1139 | 1.07 | 0.88 | 1.29 | 5.05e-01 | logit_matched_standard |
| eQTL | all_modules | 6225 | 0.2591 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4808 | 0.2552 | 0.2724 | 0.93 | 0.81 | 1.06 | 2.85e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 2314 | 0.255 | 0.2616 | 0.96 | 0.86 | 1.09 | 5.54e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 2494 | 0.2554 | 0.2616 | 0.98 | 0.87 | 1.1 | 7.41e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
