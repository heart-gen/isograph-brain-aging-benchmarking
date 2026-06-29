# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Nucleus_accumbens_basal_ganglia xQTL (gtex-aging/nucleus_accumbens_basal_ganglia)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region nucleus_accumbens_basal_ganglia`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4731 | 0.2228 | 0.2187 | 0.9 | 0.82 | 0.98 | 1.79e-02 | logit_matched |
| sQTL | pheno_sig_modules | 128 | 0.1875 | 0.2204 | 1.17 | 0.74 | 1.85 | 5.00e-01 | logit_matched |
| sQTL | go_visible_modules | 128 | 0.1875 | 0.2204 | 1.17 | 0.74 | 1.85 | 5.00e-01 | logit_matched |
| eQTL | all_modules | 5690 | 0.49 | 0.52 | 0.85 | 0.8 | 0.9 | 5.41e-07 | logit_matched |
| eQTL | pheno_sig_modules | 169 | 0.4497 | 0.5113 | 0.82 | 0.61 | 1.12 | 2.13e-01 | logit_matched |
| eQTL | go_visible_modules | 169 | 0.4497 | 0.5113 | 0.82 | 0.61 | 1.12 | 2.13e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
