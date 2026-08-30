# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4233 | 0.1857 | 0.1808 | 0.85 | 0.77 | 0.94 | 1.64e-03 | logit_matched |
| sQTL | pheno_sig_modules | 686 | 0.1822 | 0.1823 | 1.04 | 0.85 | 1.28 | 6.92e-01 | logit_matched |
| sQTL | go_invisible_modules | 443 | 0.1264 | 0.1841 | 0.69 | 0.51 | 0.92 | 1.14e-02 | logit_matched |
| sQTL | go_visible_modules | 243 | 0.284 | 0.1805 | 1.78 | 1.33 | 2.4 | 1.15e-04 | logit_matched |
| eQTL | all_modules | 5183 | 0.3556 | 0.3985 | 0.79 | 0.74 | 0.85 | 3.84e-11 | logit_matched |
| eQTL | pheno_sig_modules | 907 | 0.3374 | 0.3892 | 0.79 | 0.69 | 0.92 | 1.44e-03 | logit_matched |
| eQTL | go_invisible_modules | 610 | 0.3213 | 0.3889 | 0.73 | 0.62 | 0.87 | 4.60e-04 | logit_matched |
| eQTL | go_visible_modules | 297 | 0.3704 | 0.387 | 0.94 | 0.74 | 1.2 | 6.23e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
