# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13291 | 0.1697 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 10423 | 0.1659 | 0.1834 | 0.81 | 0.72 | 0.91 | 2.71e-04 | logit_matched |
| sQTL | go_invisible_modules | 1744 | 0.1927 | 0.1662 | 0.97 | 0.84 | 1.1 | 6.12e-01 | logit_matched |
| sQTL | go_visible_modules | 10086 | 0.1648 | 0.185 | 0.82 | 0.73 | 0.91 | 2.47e-04 | logit_matched |
| eQTL | all_modules | 17572 | 0.4033 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 13566 | 0.4017 | 0.4086 | 0.97 | 0.91 | 1.05 | 4.60e-01 | logit_matched |
| eQTL | go_invisible_modules | 2136 | 0.4167 | 0.4015 | 1.03 | 0.94 | 1.13 | 5.66e-01 | logit_matched |
| eQTL | go_visible_modules | 12991 | 0.3975 | 0.4198 | 0.92 | 0.86 | 0.98 | 1.50e-02 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
