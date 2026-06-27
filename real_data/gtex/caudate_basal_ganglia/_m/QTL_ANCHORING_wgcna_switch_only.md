# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (gtex-aging/caudate_basal_ganglia)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region caudate_basal_ganglia --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5985 | 0.212 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 134 | 0.097 | 0.2147 | 0.44 | 0.24 | 0.79 | 6.25e-03 | logit_matched |
| sQTL | go_visible_modules | 134 | 0.097 | 0.2147 | 0.44 | 0.24 | 0.79 | 6.25e-03 | logit_matched |
| eQTL | all_modules | 7229 | 0.5057 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 233 | 0.3391 | 0.5113 | 0.48 | 0.36 | 0.63 | 1.39e-07 | logit_matched |
| eQTL | go_visible_modules | 233 | 0.3391 | 0.5113 | 0.48 | 0.36 | 0.63 | 1.39e-07 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
