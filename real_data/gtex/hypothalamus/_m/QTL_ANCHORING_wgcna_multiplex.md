# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13466 | 0.1823 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 504 | 0.1627 | 0.1831 | 0.8 | 0.62 | 1.03 | 8.57e-02 | logit_matched |
| sQTL | go_invisible_modules | 151 | 0.2318 | 0.1817 | 0.89 | 0.6 | 1.33 | 5.81e-01 | logit_matched |
| sQTL | go_visible_modules | 353 | 0.1331 | 0.1836 | 0.75 | 0.55 | 1.04 | 8.51e-02 | logit_matched |
| eQTL | all_modules | 17691 | 0.3881 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 715 | 0.3259 | 0.3907 | 0.71 | 0.6 | 0.83 | 2.12e-05 | logit_matched |
| eQTL | go_invisible_modules | 158 | 0.3418 | 0.3885 | 0.71 | 0.51 | 0.99 | 4.45e-02 | logit_matched |
| eQTL | go_visible_modules | 557 | 0.3214 | 0.3902 | 0.71 | 0.59 | 0.85 | 2.06e-04 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
