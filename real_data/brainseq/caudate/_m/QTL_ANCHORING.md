# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-aging/caudate)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region caudate`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4157 | 0.2187 | 0.2301 | 0.93 | 0.84 | 1.02 | 1.08e-01 | logit_matched |
| sQTL | pheno_sig_modules | 611 | 0.234 | 0.2261 | 1.07 | 0.88 | 1.3 | 5.17e-01 | logit_matched |
| sQTL | go_invisible_modules | 611 | 0.234 | 0.2261 | 1.07 | 0.88 | 1.3 | 5.17e-01 | logit_matched |
| eQTL | all_modules | 5279 | 0.4957 | 0.536 | 0.85 | 0.8 | 0.91 | 2.02e-06 | logit_matched |
| eQTL | pheno_sig_modules | 1071 | 0.4846 | 0.5263 | 0.91 | 0.8 | 1.03 | 1.47e-01 | logit_matched |
| eQTL | go_invisible_modules | 1071 | 0.4846 | 0.5263 | 0.91 | 0.8 | 1.03 | 1.47e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
