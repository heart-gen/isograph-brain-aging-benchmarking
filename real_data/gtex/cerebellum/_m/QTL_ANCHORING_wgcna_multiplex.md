# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 9878 | 0.2671 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7565 | 0.2619 | 0.284 | 0.86 | 0.76 | 0.96 | 7.78e-03 | logit_matched |
| sQTL | go_invisible_modules | 301 | 0.3621 | 0.2641 | 1.45 | 1.12 | 1.88 | 5.33e-03 | logit_matched |
| sQTL | go_visible_modules | 7304 | 0.259 | 0.2898 | 0.83 | 0.74 | 0.92 | 6.86e-04 | logit_matched |
| eQTL | all_modules | 12647 | 0.5774 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 9186 | 0.5837 | 0.5605 | 1.09 | 1.01 | 1.18 | 2.84e-02 | logit_matched |
| eQTL | go_invisible_modules | 384 | 0.7318 | 0.5725 | 2.02 | 1.6 | 2.54 | 2.01e-09 | logit_matched |
| eQTL | go_visible_modules | 8843 | 0.578 | 0.576 | 1.0 | 0.93 | 1.08 | 9.49e-01 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
