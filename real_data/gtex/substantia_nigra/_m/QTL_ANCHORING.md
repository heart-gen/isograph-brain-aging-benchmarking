# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3992 | 0.127 | 0.1233 | 0.83 | 0.74 | 0.94 | 2.49e-03 | logit_matched |
| sQTL | pheno_sig_modules | 14 | 0.0 | 0.1245 | 0.0 | 0.0 | inf | 9.99e-01 | logit_matched |
| sQTL | go_visible_modules | 14 | 0.0 | 0.1245 | 0.0 | 0.0 | inf | 9.99e-01 | logit_matched |
| eQTL | all_modules | 4875 | 0.2556 | 0.2905 | 0.8 | 0.74 | 0.87 | 1.43e-08 | logit_matched |
| eQTL | pheno_sig_modules | 25 | 0.24 | 0.2812 | 0.79 | 0.32 | 1.99 | 6.22e-01 | logit_matched |
| eQTL | go_visible_modules | 25 | 0.24 | 0.2812 | 0.79 | 0.32 | 1.99 | 6.22e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
