# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5335 | 0.1693 | 0.1655 | 0.85 | 0.77 | 0.94 | 1.41e-03 | logit_matched |
| sQTL | pheno_sig_modules | 2298 | 0.1658 | 0.1672 | 0.87 | 0.77 | 0.99 | 3.17e-02 | logit_matched |
| sQTL | go_invisible_modules | 945 | 0.1937 | 0.165 | 1.0 | 0.84 | 1.19 | 9.98e-01 | logit_matched |
| sQTL | go_visible_modules | 1353 | 0.1463 | 0.1693 | 0.8 | 0.67 | 0.94 | 6.61e-03 | logit_matched |
| eQTL | all_modules | 6742 | 0.3558 | 0.4089 | 0.76 | 0.71 | 0.81 | 1.21e-17 | logit_matched |
| eQTL | pheno_sig_modules | 3028 | 0.3322 | 0.4007 | 0.72 | 0.66 | 0.79 | 2.02e-14 | logit_matched |
| eQTL | go_invisible_modules | 1094 | 0.372 | 0.3906 | 0.88 | 0.78 | 1.0 | 5.39e-02 | logit_matched |
| eQTL | go_visible_modules | 1934 | 0.3097 | 0.3988 | 0.67 | 0.6 | 0.74 | 1.12e-14 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
