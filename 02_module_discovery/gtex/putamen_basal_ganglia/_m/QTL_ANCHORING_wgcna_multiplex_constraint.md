# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Putamen_basal_ganglia xQTL (gtex-aging/putamen_basal_ganglia)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL], **plus gnomAD v4.1 LOEUF, missense z, and log mean expression in this cohort/region's own count matrix**. The constraint set tests whether the splicing-specificity contrast survives adjustment for selective constraint and expression level — the leading alternative explanation for module genes being cis-QTL depleted in BOTH modalities.

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region putamen_basal_ganglia --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 10759 | 0.1657 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 791 | 0.2035 | 0.1627 | 1.53 | 1.26 | 1.85 | 1.26e-05 | logit_matched_standard |
| sQTL | go_invisible_modules | 50 | 0.18 | 0.1657 | 1.27 | 0.6 | 2.7 | 5.26e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 741 | 0.2051 | 0.1628 | 1.54 | 1.27 | 1.87 | 1.50e-05 | logit_matched_standard |
| sQTL | all_modules | 10759 | 0.1657 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 791 | 0.2035 | 0.1627 | 1.35 | 1.11 | 1.64 | 2.71e-03 | logit_matched_constraint |
| sQTL | go_invisible_modules | 50 | 0.18 | 0.1657 | 1.26 | 0.58 | 2.74 | 5.55e-01 | logit_matched_constraint |
| sQTL | go_visible_modules | 741 | 0.2051 | 0.1628 | 1.35 | 1.1 | 1.65 | 3.38e-03 | logit_matched_constraint |
| eQTL | all_modules | 12909 | 0.4341 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 980 | 0.4531 | 0.4326 | 1.1 | 0.96 | 1.25 | 1.69e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 69 | 0.2609 | 0.435 | 0.46 | 0.27 | 0.79 | 5.12e-03 | logit_matched_standard |
| eQTL | go_visible_modules | 911 | 0.4676 | 0.4316 | 1.17 | 1.02 | 1.34 | 2.50e-02 | logit_matched_standard |
| eQTL | all_modules | 12909 | 0.4341 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 980 | 0.4531 | 0.4326 | 1.03 | 0.9 | 1.18 | 6.67e-01 | logit_matched_constraint |
| eQTL | go_invisible_modules | 69 | 0.2609 | 0.435 | 0.43 | 0.25 | 0.75 | 2.62e-03 | logit_matched_constraint |
| eQTL | go_visible_modules | 911 | 0.4676 | 0.4316 | 1.1 | 0.96 | 1.26 | 1.83e-01 | logit_matched_constraint |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
