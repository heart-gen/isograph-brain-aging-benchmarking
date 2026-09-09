# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3639 | 0.2399 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3023 | 0.2435 | 0.2224 | 0.94 | 0.75 | 1.17 | 5.59e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1134 | 0.2804 | 0.2216 | 1.09 | 0.92 | 1.3 | 3.14e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1889 | 0.2213 | 0.26 | 0.89 | 0.76 | 1.05 | 1.69e-01 | logit_matched_standard |
| eQTL | all_modules | 4324 | 0.5194 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3501 | 0.5301 | 0.4739 | 1.24 | 1.06 | 1.45 | 5.84e-03 | logit_matched_standard |
| eQTL | go_invisible_modules | 1277 | 0.5591 | 0.5028 | 1.23 | 1.08 | 1.41 | 2.06e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 2224 | 0.5135 | 0.5257 | 0.96 | 0.85 | 1.09 | 5.26e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
