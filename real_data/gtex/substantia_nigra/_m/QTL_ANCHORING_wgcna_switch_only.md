# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5506 | 0.1095 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4825 | 0.1086 | 0.116 | 0.81 | 0.62 | 1.05 | 1.16e-01 | logit_matched |
| sQTL | go_invisible_modules | 148 | 0.0676 | 0.1107 | 0.6 | 0.31 | 1.15 | 1.25e-01 | logit_matched |
| sQTL | go_visible_modules | 4677 | 0.1099 | 0.1074 | 0.91 | 0.71 | 1.17 | 4.70e-01 | logit_matched |
| eQTL | all_modules | 6721 | 0.257 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5778 | 0.2484 | 0.3097 | 0.75 | 0.64 | 0.87 | 1.50e-04 | logit_matched |
| eQTL | go_invisible_modules | 158 | 0.1962 | 0.2584 | 0.7 | 0.47 | 1.05 | 8.30e-02 | logit_matched |
| eQTL | go_visible_modules | 5620 | 0.2498 | 0.2934 | 0.81 | 0.7 | 0.94 | 4.32e-03 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
