# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2341 | 0.217 | 0.2257 | 0.86 | 0.77 | 0.97 | 1.06e-02 | logit_matched |
| sQTL | pheno_sig_modules | 207 | 0.2512 | 0.2238 | 1.06 | 0.76 | 1.48 | 7.13e-01 | logit_matched |
| sQTL | go_invisible_modules | 136 | 0.2868 | 0.2236 | 1.23 | 0.83 | 1.82 | 3.04e-01 | logit_matched |
| sQTL | go_visible_modules | 71 | 0.1831 | 0.2245 | 0.77 | 0.41 | 1.43 | 4.09e-01 | logit_matched |
| eQTL | all_modules | 3350 | 0.463 | 0.5326 | 0.75 | 0.7 | 0.81 | 1.57e-13 | logit_matched |
| eQTL | pheno_sig_modules | 285 | 0.407 | 0.5222 | 0.58 | 0.46 | 0.74 | 1.10e-05 | logit_matched |
| eQTL | go_invisible_modules | 165 | 0.5212 | 0.5205 | 0.91 | 0.67 | 1.24 | 5.65e-01 | logit_matched |
| eQTL | go_visible_modules | 120 | 0.25 | 0.5222 | 0.29 | 0.19 | 0.44 | 7.40e-09 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
