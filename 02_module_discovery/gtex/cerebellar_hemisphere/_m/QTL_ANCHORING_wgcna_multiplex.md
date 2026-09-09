# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cerebellar_Hemisphere xQTL (gtex-aging/cerebellar_hemisphere)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cerebellar_hemisphere --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 11890 | 0.2923 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 130 | 0.1385 | 0.294 | 0.61 | 0.36 | 1.04 | 6.72e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 130 | 0.1385 | 0.294 | 0.61 | 0.36 | 1.04 | 6.72e-02 | logit_matched_standard |
| eQTL | all_modules | 15657 | 0.6016 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 214 | 0.4065 | 0.6044 | 0.43 | 0.33 | 0.57 | 1.89e-09 | logit_matched_standard |
| eQTL | go_visible_modules | 214 | 0.4065 | 0.6044 | 0.43 | 0.33 | 0.57 | 1.89e-09 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
