# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2144 | 0.2215 | 0.2248 | 0.91 | 0.81 | 1.02 | 1.12e-01 | logit_matched |
| sQTL | pheno_sig_modules | 188 | 0.2553 | 0.2238 | 1.26 | 0.89 | 1.79 | 1.93e-01 | logit_matched |
| sQTL | go_invisible_modules | 141 | 0.2695 | 0.2238 | 1.37 | 0.92 | 2.03 | 1.23e-01 | logit_matched |
| sQTL | go_visible_modules | 47 | 0.2128 | 0.2243 | 0.97 | 0.47 | 2.01 | 9.42e-01 | logit_matched |
| eQTL | all_modules | 3074 | 0.4704 | 0.53 | 0.79 | 0.73 | 0.85 | 1.98e-09 | logit_matched |
| eQTL | pheno_sig_modules | 268 | 0.4664 | 0.5212 | 0.78 | 0.61 | 0.99 | 4.33e-02 | logit_matched |
| eQTL | go_invisible_modules | 185 | 0.5514 | 0.5202 | 1.11 | 0.83 | 1.49 | 4.89e-01 | logit_matched |
| eQTL | go_visible_modules | 83 | 0.2771 | 0.5215 | 0.33 | 0.21 | 0.54 | 9.27e-06 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
