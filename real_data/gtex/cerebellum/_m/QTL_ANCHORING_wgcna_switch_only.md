# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 804 | 0.3271 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 583 | 0.3482 | 0.2715 | 0.94 | 0.64 | 1.37 | 7.31e-01 | logit_matched |
| sQTL | go_invisible_modules | 266 | 0.3647 | 0.3086 | 0.97 | 0.68 | 1.36 | 8.39e-01 | logit_matched |
| sQTL | go_visible_modules | 317 | 0.3344 | 0.3224 | 0.98 | 0.71 | 1.37 | 9.20e-01 | logit_matched |
| eQTL | all_modules | 949 | 0.5711 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 663 | 0.5867 | 0.535 | 1.18 | 0.89 | 1.57 | 2.47e-01 | logit_matched |
| eQTL | go_invisible_modules | 311 | 0.5788 | 0.5674 | 1.03 | 0.78 | 1.35 | 8.57e-01 | logit_matched |
| eQTL | go_visible_modules | 352 | 0.5938 | 0.5578 | 1.13 | 0.87 | 1.48 | 3.61e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
