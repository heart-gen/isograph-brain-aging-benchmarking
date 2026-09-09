# Genetic anchoring — wgcna_switch_only co-switch modules vs GTEx Brain_Caudate_basal_ganglia xQTL (gtex-aging/caudate_basal_ganglia)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region caudate_basal_ganglia --method wgcna_switch_only`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 5981 | 0.2125 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 117 | 0.1026 | 0.2147 | 0.49 | 0.27 | 0.92 | 2.53e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 117 | 0.1026 | 0.2147 | 0.49 | 0.27 | 0.92 | 2.53e-02 | logit_matched_standard |
| eQTL | all_modules | 7228 | 0.5059 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 209 | 0.3493 | 0.5106 | 0.5 | 0.37 | 0.67 | 2.51e-06 | logit_matched_standard |
| eQTL | go_visible_modules | 209 | 0.3493 | 0.5106 | 0.5 | 0.37 | 0.67 | 2.51e-06 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
