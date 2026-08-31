# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3886 | 0.1225 | 0.1252 | 0.78 | 0.69 | 0.88 | 6.03e-05 | logit_matched |
| sQTL | pheno_sig_modules | 58 | 0.0862 | 0.1246 | 0.6 | 0.23 | 1.53 | 2.82e-01 | logit_matched |
| sQTL | go_invisible_modules | 30 | 0.1 | 0.1245 | 0.63 | 0.19 | 2.15 | 4.63e-01 | logit_matched |
| sQTL | go_visible_modules | 28 | 0.0714 | 0.1245 | 0.55 | 0.12 | 2.43 | 4.28e-01 | logit_matched |
| eQTL | all_modules | 4751 | 0.2572 | 0.2896 | 0.81 | 0.75 | 0.88 | 1.27e-07 | logit_matched |
| eQTL | pheno_sig_modules | 98 | 0.1633 | 0.2818 | 0.48 | 0.28 | 0.82 | 7.75e-03 | logit_matched |
| eQTL | go_invisible_modules | 40 | 0.1 | 0.2815 | 0.27 | 0.1 | 0.77 | 1.40e-02 | logit_matched |
| eQTL | go_visible_modules | 58 | 0.2069 | 0.2814 | 0.65 | 0.34 | 1.23 | 1.83e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
