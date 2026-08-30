# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6489 | 0.1581 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4412 | 0.1632 | 0.1473 | 1.0 | 0.86 | 1.16 | 9.83e-01 | logit_matched |
| sQTL | go_invisible_modules | 95 | 0.2211 | 0.1572 | 1.68 | 1.0 | 2.82 | 4.84e-02 | logit_matched |
| sQTL | go_visible_modules | 4317 | 0.1619 | 0.1506 | 0.96 | 0.83 | 1.12 | 6.29e-01 | logit_matched |
| eQTL | all_modules | 7671 | 0.3847 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5189 | 0.3849 | 0.3844 | 1.0 | 0.9 | 1.1 | 9.52e-01 | logit_matched |
| eQTL | go_invisible_modules | 112 | 0.3839 | 0.3847 | 1.03 | 0.7 | 1.52 | 8.76e-01 | logit_matched |
| eQTL | go_visible_modules | 5077 | 0.3849 | 0.3843 | 1.0 | 0.9 | 1.1 | 9.22e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
