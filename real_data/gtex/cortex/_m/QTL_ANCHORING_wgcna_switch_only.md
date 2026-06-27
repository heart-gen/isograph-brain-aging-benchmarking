# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 3628 | 0.2412 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 3013 | 0.2443 | 0.226 | 0.92 | 0.74 | 1.15 | 4.57e-01 | logit_matched |
| sQTL | go_invisible_modules | 1140 | 0.2807 | 0.2231 | 1.08 | 0.91 | 1.28 | 3.68e-01 | logit_matched |
| sQTL | go_visible_modules | 1873 | 0.2221 | 0.2615 | 0.89 | 0.76 | 1.05 | 1.64e-01 | logit_matched |
| eQTL | all_modules | 4298 | 0.5209 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3488 | 0.531 | 0.4778 | 1.23 | 1.06 | 1.44 | 7.88e-03 | logit_matched |
| eQTL | go_invisible_modules | 1298 | 0.5593 | 0.5043 | 1.23 | 1.08 | 1.41 | 1.83e-03 | logit_matched |
| eQTL | go_visible_modules | 2190 | 0.5142 | 0.528 | 0.95 | 0.84 | 1.08 | 4.41e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
