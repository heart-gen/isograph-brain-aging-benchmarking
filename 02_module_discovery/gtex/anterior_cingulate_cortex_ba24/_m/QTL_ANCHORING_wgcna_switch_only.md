# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Anterior_cingulate_cortex_BA24 xQTL (gtex-aging/anterior_cingulate_cortex_ba24)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region anterior_cingulate_cortex_ba24 --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6451 | 0.1572 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 5268 | 0.1631 | 0.131 | 1.1 | 0.91 | 1.33 | 3.28e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 2307 | 0.1886 | 0.1397 | 1.26 | 1.09 | 1.45 | 1.85e-03 | logit_matched_standard |
| sQTL | go_visible_modules | 2961 | 0.1432 | 0.1691 | 0.85 | 0.74 | 0.98 | 2.25e-02 | logit_matched_standard |
| eQTL | all_modules | 7631 | 0.3853 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 6184 | 0.3858 | 0.3829 | 0.99 | 0.88 | 1.11 | 8.70e-01 | logit_matched_standard |
| eQTL | go_invisible_modules | 2573 | 0.3956 | 0.38 | 1.03 | 0.94 | 1.14 | 5.14e-01 | logit_matched_standard |
| eQTL | go_visible_modules | 3611 | 0.3788 | 0.391 | 0.97 | 0.88 | 1.06 | 4.58e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
