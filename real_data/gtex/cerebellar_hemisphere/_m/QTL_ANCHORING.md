# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Cerebellar_Hemisphere xQTL (gtex-aging/cerebellar_hemisphere)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellar_hemisphere`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4391 | 0.3061 | 0.2863 | 0.84 | 0.77 | 0.92 | 7.89e-05 | logit_matched |
| eQTL | all_modules | 5278 | 0.6004 | 0.623 | 0.88 | 0.82 | 0.94 | 1.05e-04 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
