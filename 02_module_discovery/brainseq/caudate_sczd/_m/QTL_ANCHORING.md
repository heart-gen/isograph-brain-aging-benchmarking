# Genetic anchoring — isograph co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (brainseq-sczd)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis brainseq-sczd`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5649 | 0.2261 | 0.2236 | 0.98 | 0.9 | 1.07 | 6.21e-01 | logit_matched_standard |
| sQTL | pheno_sig_modules | 1618 | 0.2831 | 0.2167 | 1.29 | 1.15 | 1.46 | 3.29e-05 | logit_matched_standard |
| sQTL | go_invisible_modules | 1297 | 0.2791 | 0.2189 | 1.21 | 1.06 | 1.38 | 5.43e-03 | logit_matched_standard |
| sQTL | go_visible_modules | 321 | 0.2991 | 0.2228 | 1.58 | 1.23 | 2.04 | 4.13e-04 | logit_matched_standard |
| eQTL | all_modules | 6983 | 0.4875 | 0.5476 | 0.78 | 0.74 | 0.83 | 1.20e-15 | logit_matched_standard |
| eQTL | pheno_sig_modules | 1809 | 0.5412 | 0.5232 | 1.03 | 0.93 | 1.14 | 5.48e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 1427 | 0.5501 | 0.5229 | 1.06 | 0.95 | 1.18 | 3.01e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 382 | 0.5079 | 0.5253 | 0.93 | 0.76 | 1.14 | 4.97e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
