# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6283 | 0.1577 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2898 | 0.1618 | 0.1542 | 1.04 | 0.91 | 1.2 | 5.53e-01 | logit_matched |
| sQTL | go_invisible_modules | 713 | 0.1823 | 0.1546 | 1.12 | 0.91 | 1.38 | 3.00e-01 | logit_matched |
| sQTL | go_visible_modules | 2185 | 0.1551 | 0.1591 | 0.99 | 0.86 | 1.15 | 9.34e-01 | logit_matched |
| eQTL | all_modules | 7736 | 0.3692 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 3626 | 0.3593 | 0.3779 | 0.92 | 0.84 | 1.01 | 9.04e-02 | logit_matched |
| eQTL | go_invisible_modules | 845 | 0.3787 | 0.368 | 1.05 | 0.9 | 1.21 | 5.44e-01 | logit_matched |
| eQTL | go_visible_modules | 2781 | 0.3535 | 0.378 | 0.9 | 0.82 | 0.99 | 3.12e-02 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
