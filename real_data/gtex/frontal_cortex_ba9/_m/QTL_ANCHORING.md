# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6830 | 0.2151 | 0.2145 | 0.91 | 0.83 | 0.99 | 2.79e-02 | logit_matched |
| sQTL | pheno_sig_modules | 4777 | 0.2137 | 0.2154 | 0.92 | 0.84 | 1.0 | 6.02e-02 | logit_matched |
| sQTL | go_invisible_modules | 2982 | 0.2257 | 0.2118 | 1.01 | 0.91 | 1.11 | 9.02e-01 | logit_matched |
| sQTL | go_visible_modules | 1795 | 0.1939 | 0.218 | 0.83 | 0.73 | 0.95 | 5.01e-03 | logit_matched |
| eQTL | all_modules | 8211 | 0.5007 | 0.5542 | 0.79 | 0.74 | 0.83 | 1.99e-15 | logit_matched |
| eQTL | pheno_sig_modules | 5709 | 0.4957 | 0.5457 | 0.8 | 0.75 | 0.86 | 1.04e-11 | logit_matched |
| eQTL | go_invisible_modules | 3570 | 0.4888 | 0.5401 | 0.79 | 0.74 | 0.85 | 7.64e-10 | logit_matched |
| eQTL | go_visible_modules | 2139 | 0.5072 | 0.5331 | 0.9 | 0.82 | 0.99 | 2.81e-02 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
