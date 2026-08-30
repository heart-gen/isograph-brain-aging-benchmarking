# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3886 | 0.1225 | 0.1252 | 0.78 | 0.69 | 0.88 | 6.03e-05 | logit_matched |
| sQTL | pheno_sig_modules | 24 | 0.0417 | 0.1246 | 0.23 | 0.03 | 1.77 | 1.60e-01 | logit_matched |
| sQTL | go_visible_modules | 24 | 0.0417 | 0.1246 | 0.23 | 0.03 | 1.77 | 1.60e-01 | logit_matched |
| eQTL | all_modules | 4751 | 0.2572 | 0.2896 | 0.81 | 0.75 | 0.88 | 1.27e-07 | logit_matched |
| eQTL | pheno_sig_modules | 24 | 0.2083 | 0.2812 | 0.6 | 0.22 | 1.62 | 3.17e-01 | logit_matched |
| eQTL | go_visible_modules | 24 | 0.2083 | 0.2812 | 0.6 | 0.22 | 1.62 | 3.17e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
