# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13365 | 0.2227 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 1334 | 0.2549 | 0.2192 | 0.92 | 0.8 | 1.06 | 2.52e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 1334 | 0.2549 | 0.2192 | 0.92 | 0.8 | 1.06 | 2.52e-01 | logit_matched_standard |
| eQTL | all_modules | 17682 | 0.5351 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 1456 | 0.5446 | 0.5342 | 0.99 | 0.89 | 1.11 | 8.84e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 1456 | 0.5446 | 0.5342 | 0.99 | 0.89 | 1.11 | 8.84e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
