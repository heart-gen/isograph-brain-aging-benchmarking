# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Hypothalamus xQTL (gtex-aging/hypothalamus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hypothalamus --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 6125 | 0.1685 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 439 | 0.1071 | 0.1732 | 0.65 | 0.47 | 0.9 | 8.42e-03 | logit_matched_standard |
| sQTL | go_invisible_modules | 199 | 0.1005 | 0.1708 | 0.56 | 0.35 | 0.9 | 1.72e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 240 | 0.1125 | 0.1708 | 0.76 | 0.5 | 1.15 | 1.94e-01 | logit_matched_standard |
| eQTL | all_modules | 7345 | 0.3679 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 577 | 0.312 | 0.3726 | 0.77 | 0.64 | 0.93 | 5.58e-03 | logit_matched_standard |
| eQTL | go_invisible_modules | 216 | 0.2778 | 0.3706 | 0.69 | 0.51 | 0.93 | 1.54e-02 | logit_matched_standard |
| eQTL | go_visible_modules | 361 | 0.3324 | 0.3697 | 0.84 | 0.67 | 1.05 | 1.18e-01 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
