# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4514 | 0.1721 | 0.1699 | 0.91 | 0.82 | 1.01 | 6.56e-02 | logit_matched |
| sQTL | pheno_sig_modules | 3009 | 0.1562 | 0.1748 | 0.82 | 0.73 | 0.92 | 6.11e-04 | logit_matched |
| sQTL | go_invisible_modules | 1390 | 0.1727 | 0.1704 | 0.87 | 0.75 | 1.02 | 8.34e-02 | logit_matched |
| sQTL | go_visible_modules | 1619 | 0.1421 | 0.1746 | 0.81 | 0.69 | 0.94 | 6.06e-03 | logit_matched |
| eQTL | all_modules | 5485 | 0.3683 | 0.4233 | 0.78 | 0.73 | 0.83 | 2.96e-13 | logit_matched |
| eQTL | pheno_sig_modules | 3641 | 0.3483 | 0.4213 | 0.73 | 0.68 | 0.79 | 1.26e-15 | logit_matched |
| eQTL | go_invisible_modules | 1579 | 0.3768 | 0.4094 | 0.86 | 0.77 | 0.96 | 6.24e-03 | logit_matched |
| eQTL | go_visible_modules | 2062 | 0.3264 | 0.4169 | 0.68 | 0.62 | 0.76 | 3.37e-14 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
