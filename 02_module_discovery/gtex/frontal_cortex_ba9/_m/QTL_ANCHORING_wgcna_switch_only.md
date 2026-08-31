# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6553 | 0.2109 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 5113 | 0.2132 | 0.2028 | 0.98 | 0.84 | 1.14 | 8.09e-01 | logit_matched |
| sQTL | go_invisible_modules | 771 | 0.2296 | 0.2084 | 1.24 | 1.03 | 1.49 | 2.49e-02 | logit_matched |
| sQTL | go_visible_modules | 4342 | 0.2103 | 0.2121 | 0.89 | 0.78 | 1.01 | 8.20e-02 | logit_matched |
| eQTL | all_modules | 7651 | 0.5109 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5877 | 0.5125 | 0.5056 | 1.04 | 0.93 | 1.15 | 5.07e-01 | logit_matched |
| eQTL | go_invisible_modules | 858 | 0.4965 | 0.5127 | 0.96 | 0.83 | 1.1 | 5.47e-01 | logit_matched |
| eQTL | go_visible_modules | 5019 | 0.5152 | 0.5027 | 1.05 | 0.95 | 1.15 | 3.22e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
