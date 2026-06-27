# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13051 | 0.1211 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 5532 | 0.1215 | 0.1209 | 1.01 | 0.91 | 1.13 | 7.93e-01 | logit_matched |
| sQTL | go_invisible_modules | 827 | 0.1197 | 0.1212 | 0.79 | 0.63 | 0.99 | 4.41e-02 | logit_matched |
| sQTL | go_visible_modules | 5142 | 0.1208 | 0.1214 | 1.04 | 0.93 | 1.17 | 4.50e-01 | logit_matched |
| eQTL | all_modules | 17757 | 0.283 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7403 | 0.2645 | 0.2963 | 0.85 | 0.8 | 0.91 | 4.22e-06 | logit_matched |
| eQTL | go_invisible_modules | 989 | 0.2538 | 0.2848 | 0.82 | 0.71 | 0.95 | 8.85e-03 | logit_matched |
| eQTL | go_visible_modules | 6933 | 0.2628 | 0.296 | 0.85 | 0.8 | 0.91 | 5.02e-06 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
