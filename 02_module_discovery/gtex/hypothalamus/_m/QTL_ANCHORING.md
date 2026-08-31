# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4233 | 0.1857 | 0.1808 | 0.85 | 0.77 | 0.94 | 1.64e-03 | logit_matched |
| sQTL | pheno_sig_modules | 482 | 0.1349 | 0.184 | 0.65 | 0.5 | 0.86 | 2.04e-03 | logit_matched |
| sQTL | go_invisible_modules | 216 | 0.1806 | 0.1823 | 0.88 | 0.61 | 1.26 | 4.88e-01 | logit_matched |
| sQTL | go_visible_modules | 266 | 0.0977 | 0.1839 | 0.48 | 0.32 | 0.72 | 4.79e-04 | logit_matched |
| eQTL | all_modules | 5183 | 0.3556 | 0.3985 | 0.79 | 0.74 | 0.85 | 3.84e-11 | logit_matched |
| eQTL | pheno_sig_modules | 608 | 0.3026 | 0.3895 | 0.64 | 0.54 | 0.77 | 9.22e-07 | logit_matched |
| eQTL | go_invisible_modules | 234 | 0.3333 | 0.3874 | 0.75 | 0.57 | 0.99 | 4.10e-02 | logit_matched |
| eQTL | go_visible_modules | 374 | 0.2834 | 0.3888 | 0.59 | 0.47 | 0.74 | 4.73e-06 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
