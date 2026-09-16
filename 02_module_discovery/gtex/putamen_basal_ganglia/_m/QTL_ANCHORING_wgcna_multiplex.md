# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Putamen_basal_ganglia xQTL (gtex-aging/putamen_basal_ganglia)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region putamen_basal_ganglia --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 12633 | 0.181 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 935 | 0.2086 | 0.1788 | 1.4 | 1.18 | 1.66 | 1.42e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 52 | 0.1731 | 0.1811 | 1.11 | 0.53 | 2.32 | 7.90e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 883 | 0.2106 | 0.1788 | 1.41 | 1.18 | 1.68 | 1.24e-04 | logit_matched_standard |
| eQTL | all_modules | 16543 | 0.4577 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 1191 | 0.4576 | 0.4577 | 1.01 | 0.9 | 1.14 | 8.31e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 75 | 0.28 | 0.4585 | 0.46 | 0.28 | 0.77 | 2.85e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 1116 | 0.4695 | 0.4568 | 1.07 | 0.94 | 1.21 | 2.97e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
