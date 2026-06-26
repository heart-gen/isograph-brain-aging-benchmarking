# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 4242 | 0.2254 | 0.2101 | 0.96 | 0.88 | 1.06 | 4.12e-01 | logit_matched |
| sQTL | pheno_sig_modules | 3042 | 0.2133 | 0.2152 | 0.86 | 0.77 | 0.95 | 4.22e-03 | logit_matched |
| sQTL | go_invisible_modules | 1778 | 0.248 | 0.2098 | 1.04 | 0.92 | 1.17 | 5.49e-01 | logit_matched |
| sQTL | go_visible_modules | 1264 | 0.1646 | 0.2199 | 0.67 | 0.57 | 0.79 | 1.30e-06 | logit_matched |
| eQTL | all_modules | 5218 | 0.4994 | 0.5424 | 0.82 | 0.77 | 0.88 | 3.05e-09 | logit_matched |
| eQTL | pheno_sig_modules | 3752 | 0.4915 | 0.5401 | 0.81 | 0.75 | 0.87 | 5.77e-09 | logit_matched |
| eQTL | go_invisible_modules | 2043 | 0.5056 | 0.5332 | 0.85 | 0.78 | 0.94 | 8.47e-04 | logit_matched |
| eQTL | go_visible_modules | 1709 | 0.4745 | 0.5358 | 0.8 | 0.72 | 0.88 | 8.79e-06 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
