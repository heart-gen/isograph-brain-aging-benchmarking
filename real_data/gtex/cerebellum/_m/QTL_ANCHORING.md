# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Cerebellum xQTL (gtex-aging/cerebellum)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellum`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 2173 | 0.3212 | 0.2747 | 0.96 | 0.86 | 1.07 | 4.79e-01 | logit_matched |
| sQTL | pheno_sig_modules | 637 | 0.3014 | 0.2815 | 0.81 | 0.67 | 0.98 | 3.00e-02 | logit_matched |
| sQTL | go_invisible_modules | 459 | 0.3769 | 0.279 | 0.93 | 0.76 | 1.14 | 4.77e-01 | logit_matched |
| sQTL | go_visible_modules | 178 | 0.1067 | 0.2848 | 0.44 | 0.27 | 0.72 | 1.22e-03 | logit_matched |
| eQTL | all_modules | 2578 | 0.5877 | 0.6142 | 0.85 | 0.78 | 0.92 | 1.49e-04 | logit_matched |
| eQTL | pheno_sig_modules | 847 | 0.5478 | 0.6135 | 0.71 | 0.62 | 0.82 | 2.33e-06 | logit_matched |
| eQTL | go_invisible_modules | 493 | 0.6592 | 0.6092 | 1.16 | 0.96 | 1.4 | 1.24e-01 | logit_matched |
| eQTL | go_visible_modules | 354 | 0.3927 | 0.6147 | 0.38 | 0.31 | 0.48 | 4.37e-18 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
