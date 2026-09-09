# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Frontal_Cortex_BA9 xQTL (gtex-aging/frontal_cortex_ba9)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region frontal_cortex_ba9 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6537 | 0.2082 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 5121 | 0.2121 | 0.1942 | 1.02 | 0.87 | 1.19 | 8.10e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 2400 | 0.2329 | 0.1939 | 1.15 | 1.01 | 1.31 | 3.20e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 2721 | 0.1937 | 0.2186 | 0.88 | 0.78 | 1.0 | 5.50e-02 | logit_matched_standard |
| eQTL | all_modules | 7646 | 0.5107 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5890 | 0.5115 | 0.508 | 1.01 | 0.91 | 1.13 | 8.06e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 2656 | 0.5173 | 0.5072 | 1.02 | 0.93 | 1.12 | 6.98e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 3234 | 0.5068 | 0.5136 | 0.99 | 0.91 | 1.09 | 8.70e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
