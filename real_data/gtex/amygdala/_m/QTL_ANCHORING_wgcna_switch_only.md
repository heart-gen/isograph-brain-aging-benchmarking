# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Amygdala xQTL (gtex-aging/amygdala)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region amygdala --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4858 | 0.1136 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3751 | 0.109 | 0.1292 | 0.94 | 0.76 | 1.16 | 5.45e-01 | logit_matched |
| sQTL | go_invisible_modules | 1982 | 0.1054 | 0.1193 | 0.88 | 0.73 | 1.06 | 1.82e-01 | logit_matched |
| sQTL | go_visible_modules | 1769 | 0.1131 | 0.114 | 1.09 | 0.9 | 1.31 | 4.04e-01 | logit_matched |
| eQTL | all_modules | 6069 | 0.2595 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4591 | 0.2566 | 0.2686 | 0.96 | 0.84 | 1.1 | 5.36e-01 | logit_matched |
| eQTL | go_invisible_modules | 2247 | 0.2555 | 0.2619 | 0.97 | 0.86 | 1.09 | 6.30e-01 | logit_matched |
| eQTL | go_visible_modules | 2344 | 0.2577 | 0.2607 | 1.0 | 0.88 | 1.12 | 9.47e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
