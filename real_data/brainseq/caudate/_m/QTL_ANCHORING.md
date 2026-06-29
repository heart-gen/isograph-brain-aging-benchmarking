# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-aging/caudate)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-aging --region caudate`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4254 | 0.2151 | 0.232 | 0.91 | 0.83 | 1.0 | 4.82e-02 | logit_matched |
| sQTL | pheno_sig_modules | 646 | 0.2353 | 0.226 | 1.07 | 0.88 | 1.29 | 5.22e-01 | logit_matched |
| sQTL | go_invisible_modules | 646 | 0.2353 | 0.226 | 1.07 | 0.88 | 1.29 | 5.22e-01 | logit_matched |
| eQTL | all_modules | 5465 | 0.4917 | 0.5385 | 0.84 | 0.78 | 0.89 | 4.79e-08 | logit_matched |
| eQTL | pheno_sig_modules | 1175 | 0.4774 | 0.5271 | 0.89 | 0.79 | 1.0 | 4.97e-02 | logit_matched |
| eQTL | go_invisible_modules | 1175 | 0.4774 | 0.5271 | 0.89 | 0.79 | 1.0 | 4.97e-02 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
