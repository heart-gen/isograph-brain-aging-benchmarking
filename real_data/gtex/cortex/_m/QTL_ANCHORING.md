# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2342 | 0.2605 | 0.2177 | 1.06 | 0.95 | 1.18 | 2.92e-01 | logit_matched |
| sQTL | pheno_sig_modules | 1230 | 0.2634 | 0.2212 | 1.02 | 0.88 | 1.17 | 8.10e-01 | logit_matched |
| sQTL | go_invisible_modules | 1041 | 0.2786 | 0.2206 | 1.04 | 0.89 | 1.21 | 6.27e-01 | logit_matched |
| sQTL | go_visible_modules | 189 | 0.1799 | 0.2257 | 0.89 | 0.6 | 1.32 | 5.54e-01 | logit_matched |
| eQTL | all_modules | 2921 | 0.5043 | 0.5434 | 0.83 | 0.76 | 0.9 | 2.84e-06 | logit_matched |
| eQTL | pheno_sig_modules | 1438 | 0.4993 | 0.5404 | 0.81 | 0.73 | 0.91 | 2.23e-04 | logit_matched |
| eQTL | go_invisible_modules | 1141 | 0.5092 | 0.539 | 0.82 | 0.73 | 0.93 | 1.73e-03 | logit_matched |
| eQTL | go_visible_modules | 297 | 0.4613 | 0.5384 | 0.8 | 0.63 | 1.01 | 5.81e-02 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
