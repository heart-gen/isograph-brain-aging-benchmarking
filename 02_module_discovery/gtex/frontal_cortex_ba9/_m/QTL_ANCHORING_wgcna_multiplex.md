# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13415 | 0.2141 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 8025 | 0.2204 | 0.2046 | 1.03 | 0.94 | 1.12 | 5.72e-01 | logit_matched |
| sQTL | go_invisible_modules | 1683 | 0.2436 | 0.2099 | 0.99 | 0.87 | 1.12 | 8.33e-01 | logit_matched |
| sQTL | go_visible_modules | 7320 | 0.2178 | 0.2097 | 1.05 | 0.96 | 1.14 | 3.06e-01 | logit_matched |
| eQTL | all_modules | 17569 | 0.5273 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 10046 | 0.533 | 0.5197 | 1.06 | 1.0 | 1.12 | 6.70e-02 | logit_matched |
| eQTL | go_invisible_modules | 1883 | 0.5236 | 0.5278 | 0.95 | 0.86 | 1.04 | 2.82e-01 | logit_matched |
| eQTL | go_visible_modules | 9249 | 0.5355 | 0.5183 | 1.08 | 1.02 | 1.15 | 8.10e-03 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
