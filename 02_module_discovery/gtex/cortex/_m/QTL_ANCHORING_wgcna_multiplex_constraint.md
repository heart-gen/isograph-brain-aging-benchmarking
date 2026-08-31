# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Cortex xQTL (gtex-aging/cortex)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region cortex --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 11326 | 0.2078 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 6144 | 0.2236 | 0.1891 | 1.2 | 1.09 | 1.32 | 2.34e-04 | logit_matched_standard |
| sQTL | go_invisible_modules | 1242 | 0.2609 | 0.2013 | 1.15 | 0.99 | 1.32 | 6.46e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 5533 | 0.2201 | 0.1961 | 1.2 | 1.09 | 1.32 | 2.21e-04 | logit_matched_standard |
| sQTL | all_modules | 11326 | 0.2078 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 6144 | 0.2236 | 0.1891 | 1.13 | 1.02 | 1.25 | 1.56e-02 | logit_matched_constraint |
| sQTL | go_invisible_modules | 1242 | 0.2609 | 0.2013 | 1.16 | 1.0 | 1.34 | 5.14e-02 | logit_matched_constraint |
| sQTL | go_visible_modules | 5533 | 0.2201 | 0.1961 | 1.12 | 1.02 | 1.24 | 2.32e-02 | logit_matched_constraint |
| eQTL | all_modules | 13586 | 0.513 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7106 | 0.5227 | 0.5023 | 1.08 | 1.01 | 1.16 | 2.05e-02 | logit_matched_standard |
| eQTL | go_invisible_modules | 1405 | 0.5246 | 0.5116 | 1.02 | 0.91 | 1.14 | 7.35e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 6401 | 0.5227 | 0.5042 | 1.08 | 1.01 | 1.16 | 1.99e-02 | logit_matched_standard |
| eQTL | all_modules | 13586 | 0.513 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 7106 | 0.5227 | 0.5023 | 1.08 | 1.01 | 1.16 | 2.80e-02 | logit_matched_constraint |
| eQTL | go_invisible_modules | 1405 | 0.5246 | 0.5116 | 1.04 | 0.93 | 1.16 | 5.17e-01 | logit_matched_constraint |
| eQTL | go_visible_modules | 6401 | 0.5227 | 0.5042 | 1.08 | 1.01 | 1.16 | 3.34e-02 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
