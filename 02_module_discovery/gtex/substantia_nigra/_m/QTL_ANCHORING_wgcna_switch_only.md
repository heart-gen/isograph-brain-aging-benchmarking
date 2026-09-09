# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Substantia_nigra xQTL (gtex-aging/substantia_nigra)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region substantia_nigra --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5520 | 0.1085 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4793 | 0.1079 | 0.1128 | 0.84 | 0.65 | 1.08 | 1.78e-01 | logit_matched_standard |
| sQTL | go_invisible_modules | 109 | 0.0734 | 0.1092 | 0.66 | 0.32 | 1.39 | 2.75e-01 | logit_matched_standard |
| sQTL | go_visible_modules | 4684 | 0.1087 | 0.1077 | 0.9 | 0.71 | 1.15 | 4.08e-01 | logit_matched_standard |
| eQTL | all_modules | 6742 | 0.2578 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 5746 | 0.2483 | 0.3122 | 0.74 | 0.64 | 0.85 | 5.15e-05 | logit_matched_standard |
| eQTL | go_invisible_modules | 115 | 0.1565 | 0.2595 | 0.53 | 0.32 | 0.88 | 1.44e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 5631 | 0.2502 | 0.2961 | 0.8 | 0.7 | 0.93 | 2.58e-03 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
