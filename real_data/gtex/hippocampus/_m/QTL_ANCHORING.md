# Genetic anchoring — IsoGraph co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with IsoGraph's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5574 | 0.1679 | 0.1663 | 0.85 | 0.77 | 0.94 | 1.19e-03 | logit_matched |
| sQTL | pheno_sig_modules | 1875 | 0.1563 | 0.1687 | 0.84 | 0.73 | 0.97 | 1.54e-02 | logit_matched |
| sQTL | go_invisible_modules | 590 | 0.1814 | 0.1663 | 1.02 | 0.81 | 1.27 | 8.86e-01 | logit_matched |
| sQTL | go_visible_modules | 1285 | 0.1447 | 0.1694 | 0.77 | 0.65 | 0.92 | 3.07e-03 | logit_matched |
| eQTL | all_modules | 7118 | 0.3581 | 0.4093 | 0.77 | 0.72 | 0.82 | 2.19e-16 | logit_matched |
| eQTL | pheno_sig_modules | 2563 | 0.3262 | 0.3997 | 0.72 | 0.66 | 0.78 | 2.50e-13 | logit_matched |
| eQTL | go_invisible_modules | 778 | 0.3728 | 0.3902 | 0.91 | 0.78 | 1.06 | 2.14e-01 | logit_matched |
| eQTL | go_visible_modules | 1785 | 0.3059 | 0.3984 | 0.66 | 0.59 | 0.73 | 1.51e-14 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
