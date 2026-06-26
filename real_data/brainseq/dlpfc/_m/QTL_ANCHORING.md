# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (brainseq-aging/dlpfc)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region dlpfc`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3869 | 0.2024 | 0.2244 | 0.86 | 0.78 | 0.95 | 2.66e-03 | logit_matched |
| eQTL | all_modules | 4406 | 0.4855 | 0.5264 | 0.83 | 0.77 | 0.89 | 1.21e-07 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
