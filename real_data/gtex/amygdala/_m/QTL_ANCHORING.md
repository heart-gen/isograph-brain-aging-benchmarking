# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3387 | 0.1193 | 0.1219 | 0.77 | 0.67 | 0.87 | 4.27e-05 | logit_matched |
| sQTL | pheno_sig_modules | 1168 | 0.0942 | 0.1238 | 0.61 | 0.49 | 0.75 | 4.44e-06 | logit_matched |
| sQTL | go_invisible_modules | 491 | 0.0896 | 0.1224 | 0.63 | 0.45 | 0.87 | 4.72e-03 | logit_matched |
| sQTL | go_visible_modules | 677 | 0.0975 | 0.1225 | 0.63 | 0.48 | 0.82 | 7.17e-04 | logit_matched |
| eQTL | all_modules | 4240 | 0.2592 | 0.2942 | 0.81 | 0.75 | 0.88 | 1.73e-07 | logit_matched |
| eQTL | pheno_sig_modules | 1528 | 0.2271 | 0.2915 | 0.7 | 0.62 | 0.8 | 3.65e-08 | logit_matched |
| eQTL | go_invisible_modules | 551 | 0.216 | 0.2882 | 0.67 | 0.54 | 0.82 | 1.41e-04 | logit_matched |
| eQTL | go_visible_modules | 977 | 0.2334 | 0.289 | 0.74 | 0.64 | 0.86 | 1.28e-04 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
