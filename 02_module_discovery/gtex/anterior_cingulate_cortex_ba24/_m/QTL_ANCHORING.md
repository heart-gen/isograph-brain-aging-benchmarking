# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4757 | 0.1682 | 0.172 | 0.87 | 0.79 | 0.96 | 5.24e-03 | logit_matched |
| sQTL | pheno_sig_modules | 3078 | 0.1514 | 0.1764 | 0.75 | 0.67 | 0.84 | 1.24e-06 | logit_matched |
| sQTL | go_invisible_modules | 1392 | 0.1667 | 0.1711 | 0.85 | 0.73 | 1.0 | 4.53e-02 | logit_matched |
| sQTL | go_visible_modules | 1686 | 0.1388 | 0.1752 | 0.72 | 0.62 | 0.84 | 2.06e-05 | logit_matched |
| eQTL | all_modules | 5799 | 0.3656 | 0.426 | 0.76 | 0.71 | 0.81 | 3.45e-16 | logit_matched |
| eQTL | pheno_sig_modules | 3743 | 0.354 | 0.4203 | 0.74 | 0.69 | 0.8 | 1.75e-14 | logit_matched |
| eQTL | go_invisible_modules | 1566 | 0.3563 | 0.4114 | 0.77 | 0.69 | 0.86 | 3.67e-06 | logit_matched |
| eQTL | go_visible_modules | 2177 | 0.3523 | 0.414 | 0.77 | 0.7 | 0.84 | 4.44e-08 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
