# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4284 | 0.1872 | 0.1801 | 0.88 | 0.8 | 0.97 | 1.11e-02 | logit_matched |
| sQTL | pheno_sig_modules | 731 | 0.1683 | 0.1831 | 0.89 | 0.73 | 1.09 | 2.65e-01 | logit_matched |
| sQTL | go_invisible_modules | 545 | 0.1908 | 0.1819 | 1.05 | 0.84 | 1.31 | 6.92e-01 | logit_matched |
| sQTL | go_visible_modules | 186 | 0.1022 | 0.1834 | 0.5 | 0.31 | 0.81 | 4.84e-03 | logit_matched |
| eQTL | all_modules | 5262 | 0.3601 | 0.397 | 0.82 | 0.77 | 0.88 | 9.21e-09 | logit_matched |
| eQTL | pheno_sig_modules | 907 | 0.3374 | 0.3892 | 0.78 | 0.67 | 0.9 | 4.93e-04 | logit_matched |
| eQTL | go_invisible_modules | 623 | 0.3579 | 0.3877 | 0.87 | 0.74 | 1.03 | 1.10e-01 | logit_matched |
| eQTL | go_visible_modules | 284 | 0.2923 | 0.3882 | 0.61 | 0.47 | 0.79 | 1.86e-04 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
