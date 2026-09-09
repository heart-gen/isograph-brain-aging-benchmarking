# Genetic anchoring — wgcna_multiplex co-switch modules vs GTEx Brain_Hippocampus xQTL (gtex-aging/hippocampus)

Power-matched enrichment (logistic: sGene/eGene status ~ module membership + covariates) within each xQTL's tested-gene universe intersected with the method's tested genes.

Covariates: log cis-variant count, log gene length, log isoform count [+ log intron group size for sQTL].

Reproduce: `python -m isograph_benchmark.real_data.qtl_anchoring --analysis gtex-aging --region hippocampus --method wgcna_multiplex`

## Matched odds ratios

| xqtl_kind | module_set | n_foreground | rate_fg | rate_bg | odds_ratio | or_ci_low | or_ci_high | pvalue | fit_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sQTL | all_modules | 12913 | 0.1675 | nan | nan | nan | nan | NA | insufficient |
| sQTL | pheno_sig_modules | 4743 | 0.1516 | 0.1767 | 0.88 | 0.8 | 0.98 | 1.50e-02 | logit_matched_standard |
| sQTL | go_invisible_modules | 216 | 0.1111 | 0.1685 | 0.62 | 0.4 | 0.97 | 3.53e-02 | logit_matched_standard |
| sQTL | go_visible_modules | 4591 | 0.1527 | 0.1757 | 0.9 | 0.81 | 0.99 | 3.35e-02 | logit_matched_standard |
| eQTL | all_modules | 17225 | 0.3904 | nan | nan | nan | nan | NA | insufficient |
| eQTL | pheno_sig_modules | 6491 | 0.3445 | 0.4182 | 0.73 | 0.69 | 0.78 | 5.50e-21 | logit_matched_standard |
| eQTL | go_invisible_modules | 291 | 0.2715 | 0.3925 | 0.57 | 0.44 | 0.75 | 3.24e-05 | logit_matched_standard |
| eQTL | go_visible_modules | 6268 | 0.3468 | 0.4154 | 0.75 | 0.7 | 0.8 | 4.92e-18 | logit_matched_standard |

## Reading

- **Read the two arms together, not the ratio alone.** Co-switch module genes are cis-QTL *depleted* for BOTH modalities (OR < 1 in every tissue tested); the result is that splicing-QTL is spared RELATIVE to expression-QTL, i.e. sQTL OR / eQTL OR > 1. It is a ratio of two depletions, not an enrichment, and must never be described as sQTL enrichment.
- Matching on cis-variant count / gene length / isoform multiplicity controls the dominant QTL-detectability confound; `rate_fg` vs `rate_bg` is the raw (unmatched) contrast for reference.
- Scope: cis-sQTL enrichment shows module *members* undergo genetically regulated splicing; it does not by itself prove the *co-switching* is genetic (a shared trans regulator / cell composition could coordinate it).
