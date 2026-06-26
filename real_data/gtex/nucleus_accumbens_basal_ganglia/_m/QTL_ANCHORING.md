# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Nucleus_accumbens_basal_ganglia xQTL (gtex-aging/nucleus_accumbens_basal_ganglia)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region nucleus_accumbens_basal_ganglia`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4357 | 0.2203 | 0.22 | 0.86 | 0.79 | 0.95 | 2.11e-03 | logit_matched |
| sQTL | pheno_sig_modules | 80 | 0.0375 | 0.2212 | 0.14 | 0.04 | 0.45 | 9.38e-04 | logit_matched |
| sQTL | go_visible_modules | 80 | 0.0375 | 0.2212 | 0.14 | 0.04 | 0.45 | 9.38e-04 | logit_matched |
| eQTL | all_modules | 5261 | 0.4942 | 0.5174 | 0.87 | 0.82 | 0.93 | 5.32e-05 | logit_matched |
| eQTL | pheno_sig_modules | 150 | 0.32 | 0.5123 | 0.42 | 0.3 | 0.6 | 1.03e-06 | logit_matched |
| eQTL | go_visible_modules | 150 | 0.32 | 0.5123 | 0.42 | 0.3 | 0.6 | 1.03e-06 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
