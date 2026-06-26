# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Putamen_basal_ganglia xQTL (gtex-aging/putamen_basal_ganglia)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region putamen_basal_ganglia`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4199 | 0.1836 | 0.1793 | 0.84 | 0.76 | 0.93 | 1.01e-03 | logit_matched |
| sQTL | pheno_sig_modules | 282 | 0.1667 | 0.181 | 0.86 | 0.62 | 1.19 | 3.55e-01 | logit_matched |
| sQTL | go_invisible_modules | 151 | 0.2119 | 0.1803 | 1.13 | 0.75 | 1.69 | 5.56e-01 | logit_matched |
| sQTL | go_visible_modules | 131 | 0.1145 | 0.1813 | 0.57 | 0.33 | 0.99 | 4.51e-02 | logit_matched |
| eQTL | all_modules | 5172 | 0.4397 | 0.4682 | 0.84 | 0.78 | 0.9 | 2.78e-07 | logit_matched |
| eQTL | pheno_sig_modules | 408 | 0.3701 | 0.4621 | 0.65 | 0.53 | 0.79 | 2.92e-05 | logit_matched |
| eQTL | go_invisible_modules | 182 | 0.4451 | 0.4601 | 0.93 | 0.69 | 1.24 | 6.05e-01 | logit_matched |
| eQTL | go_visible_modules | 226 | 0.3097 | 0.4619 | 0.48 | 0.36 | 0.64 | 4.32e-07 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
