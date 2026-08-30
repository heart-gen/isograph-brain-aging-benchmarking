# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Putamen_basal_ganglia xQTL (gtex-aging/putamen_basal_ganglia)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region putamen_basal_ganglia`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5354 | 0.1802 | 0.1809 | 0.85 | 0.77 | 0.94 | 9.67e-04 | logit_matched |
| sQTL | pheno_sig_modules | 548 | 0.1843 | 0.1805 | 1.01 | 0.8 | 1.26 | 9.61e-01 | logit_matched |
| sQTL | go_invisible_modules | 216 | 0.2361 | 0.1797 | 1.31 | 0.94 | 1.82 | 1.10e-01 | logit_matched |
| sQTL | go_visible_modules | 332 | 0.1506 | 0.1814 | 0.81 | 0.59 | 1.11 | 1.98e-01 | logit_matched |
| eQTL | all_modules | 6423 | 0.43 | 0.4767 | 0.78 | 0.74 | 0.83 | 3.78e-14 | logit_matched |
| eQTL | pheno_sig_modules | 656 | 0.4085 | 0.4619 | 0.82 | 0.7 | 0.96 | 1.22e-02 | logit_matched |
| eQTL | go_invisible_modules | 250 | 0.488 | 0.4596 | 1.11 | 0.86 | 1.42 | 4.22e-01 | logit_matched |
| eQTL | go_visible_modules | 406 | 0.3596 | 0.4623 | 0.67 | 0.55 | 0.83 | 1.54e-04 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
