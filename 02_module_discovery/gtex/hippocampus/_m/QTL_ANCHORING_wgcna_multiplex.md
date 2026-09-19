# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 12773 | 0.1666 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 2792 | 0.1357 | 0.1752 | 0.79 | 0.7 | 0.9 | 2.28e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 267 | 0.1723 | 0.1665 | 0.92 | 0.66 | 1.28 | 6.21e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 2601 | 0.1319 | 0.1755 | 0.77 | 0.67 | 0.87 | 5.32e-05 | logit_matched_standard |
| eQTL | all_modules | 16698 | 0.3894 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 4002 | 0.3143 | 0.413 | 0.66 | 0.61 | 0.71 | 3.65e-27 | logit_matched_standard |
| eQTL | go_invisible_modules | 339 | 0.3097 | 0.391 | 0.7 | 0.55 | 0.88 | 2.48e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 3742 | 0.3135 | 0.4113 | 0.66 | 0.61 | 0.71 | 1.28e-25 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
