# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Nucleus_accumbens_basal_ganglia xQTL (gtex-aging/nucleus_accumbens_basal_ganglia)

Power-matched enrichment (logistic: qtl status ~ module membership + log cis-variant count + log gene length + log isoform count [+ log intron group size for sQTL]) within each xQTL's tested-gene universe intersected with the method's tested genes.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region nucleus_accumbens_basal_ganglia --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 12981 | 0.2193 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 225 | 0.1422 | 0.2207 | 0.67 | 0.45 | 0.99 | 4.53e-02 | logit_matched |
| sQTL | go_visible_modules | 225 | 0.1422 | 0.2207 | 0.67 | 0.45 | 0.99 | 4.53e-02 | logit_matched |
| eQTL | all_modules | 17088 | 0.5088 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 420 | 0.4 | 0.5115 | 0.61 | 0.5 | 0.75 | 1.25e-06 | logit_matched |
| eQTL | go_visible_modules | 420 | 0.4 | 0.5115 | 0.61 | 0.5 | 0.75 | 1.25e-06 | logit_matched |

## Reading

- On-thesis result: **sQTL OR > 1 and significant** while the matched **eQTL OR is near 1** for the same module set => the genetic signal on co-switch modules is splicing-specific (DTU-without-DGE) rather than expression-level.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
