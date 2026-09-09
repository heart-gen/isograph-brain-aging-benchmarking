# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 13397 | 0.214 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 7996 | 0.2209 | 0.2039 | 1.03 | 0.94 | 1.12 | 5.49e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 1682 | 0.2432 | 0.2098 | 0.98 | 0.86 | 1.11 | 7.23e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 7275 | 0.2179 | 0.2094 | 1.05 | 0.96 | 1.14 | 3.20e-01 | logit_matched_standard |
| eQTL | all_modules | 17561 | 0.5268 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 10011 | 0.5325 | 0.5193 | 1.06 | 0.99 | 1.12 | 7.27e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 1873 | 0.528 | 0.5267 | 0.97 | 0.88 | 1.06 | 4.89e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 9201 | 0.5343 | 0.5187 | 1.08 | 1.01 | 1.14 | 1.47e-02 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
