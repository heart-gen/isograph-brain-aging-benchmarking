# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 12905 | 0.1676 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4267 | 0.1362 | 0.1831 | 0.76 | 0.69 | 0.85 | 5.83e-07 | logit_matched |
| sQTL | go_invisible_modules | 217 | 0.1152 | 0.1685 | 0.65 | 0.42 | 1.0 | 5.05e-02 | logit_matched |
| sQTL | go_visible_modules | 4117 | 0.1368 | 0.1821 | 0.77 | 0.69 | 0.86 | 2.15e-06 | logit_matched |
| eQTL | all_modules | 17207 | 0.3904 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 6012 | 0.3395 | 0.4178 | 0.73 | 0.68 | 0.77 | 1.05e-21 | logit_matched |
| eQTL | go_invisible_modules | 287 | 0.2683 | 0.3925 | 0.56 | 0.43 | 0.73 | 2.09e-05 | logit_matched |
| eQTL | go_visible_modules | 5796 | 0.3421 | 0.415 | 0.74 | 0.7 | 0.79 | 1.58e-18 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
