# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3065 | 0.1171 | 0.1224 | 0.74 | 0.64 | 0.84 | 6.36e-06 | logit_matched |
| sQTL | pheno_sig_modules | 1159 | 0.0984 | 0.1234 | 0.66 | 0.53 | 0.81 | 9.91e-05 | logit_matched |
| sQTL | go_invisible_modules | 515 | 0.0971 | 0.1222 | 0.71 | 0.52 | 0.96 | 2.71e-02 | logit_matched |
| sQTL | go_visible_modules | 644 | 0.0994 | 0.1223 | 0.65 | 0.5 | 0.86 | 2.32e-03 | logit_matched |
| eQTL | all_modules | 3840 | 0.2531 | 0.2949 | 0.78 | 0.72 | 0.85 | 3.26e-09 | logit_matched |
| eQTL | pheno_sig_modules | 1507 | 0.2236 | 0.2917 | 0.69 | 0.61 | 0.79 | 1.47e-08 | logit_matched |
| eQTL | go_invisible_modules | 575 | 0.2243 | 0.2881 | 0.71 | 0.58 | 0.87 | 7.87e-04 | logit_matched |
| eQTL | go_visible_modules | 932 | 0.2232 | 0.2894 | 0.7 | 0.6 | 0.82 | 1.23e-05 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
