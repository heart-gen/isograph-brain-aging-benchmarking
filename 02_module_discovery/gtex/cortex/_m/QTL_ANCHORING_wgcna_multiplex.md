# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13376 | 0.2234 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7402 | 0.2348 | 0.2092 | 1.15 | 1.06 | 1.26 | 1.36e-03 | logit_matched_standard |
| sQTL | go_invisible_modules | 1413 | 0.2767 | 0.2171 | 1.09 | 0.96 | 1.24 | 2.01e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 6698 | 0.2308 | 0.2159 | 1.16 | 1.07 | 1.26 | 6.50e-04 | logit_matched_standard |
| eQTL | all_modules | 17851 | 0.5338 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 9315 | 0.5487 | 0.5175 | 1.13 | 1.07 | 1.2 | 3.46e-05 | logit_matched_standard |
| eQTL | go_invisible_modules | 1684 | 0.5487 | 0.5322 | 1.02 | 0.92 | 1.13 | 6.68e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 8443 | 0.5491 | 0.52 | 1.13 | 1.07 | 1.2 | 2.85e-05 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
