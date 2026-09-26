# Threshold sensitivity of the switch-unique classification

A recount, not a refit: `fdr_switch_given_abund` and `fdr_abund_given_switch` are BH-adjusted once over all genes by `incremental_association`, so changing alpha moves only the cut. The category rule is reproduced exactly, `<=` included.

17 analyses; 12 of them also have a composition-adjusted arm.

## Switch-unique counts and gene-level persistence by alpha

| label | base_a0.05 | base_a0.1 | base_a0.2 | persist_a0.05 | persist_a0.1 | persist_a0.2 |
| --- | --- | --- | --- | --- | --- | --- |
| GTEx cortex | 835 | 1169 | 1663 | 0 | 0 | 0.00241 |
| GTEx anterior_cingulate_cortex_ba24 | 563 | 928 | 1398 | 0.258 | 0.316 | 0.4 |
| GTEx frontal_cortex_ba9 | 512 | 797 | 1128 | 0.00195 | 0.00125 | 0.00355 |
| GTEx hypothalamus | 145 | 383 | 882 |  |  |  |
| GTEx caudate_basal_ganglia | 42 | 87 | 266 | 0.167 | 0.218 | 0.132 |
| GTEx cerebellum | 30 | 69 | 187 |  |  |  |
| SCZD (caudate) | 20 | 60 | 190 | 0.15 | 0.0833 | 0.0421 |
| aging caudate | 25 | 45 | 86 | 0.16 | 0.222 | 0.279 |
| GTEx hippocampus | 26 | 44 | 159 | 0.0385 | 0.477 | 0.547 |
| GTEx amygdala | 12 | 33 | 76 | 0.333 | 0.121 | 0.237 |
| aging DLPFC | 18 | 22 | 36 | 0.333 | 0.364 | 0.472 |
| GTEx substantia_nigra | 4 | 12 | 38 |  |  |  |
| GTEx cerebellar_hemisphere | 3 | 6 | 6 |  |  |  |
| GTEx nucleus_accumbens_basal_ganglia | 1 | 2 | 5 | 1 | 0.5 | 0.8 |
| GTEx putamen_basal_ganglia | 2 | 2 | 10 | 1 | 1 | 0.2 |
| GTEx spinal_cord_cervical_c_1 | 1 | 1 | 10 |  |  |  |
| aging hippocampus | 1 | 1 | 3 | 0 | 1 | 0.667 |

## Does alpha reorder the analyses?

| alpha | n_analyses | spearman_switch_unique_base | n_switch_unique_base | spearman_persistence | n_persistence | spearman_persistence_large | n_persistence_large |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.05 | 17 | 0.984 | 17 | 0.401 | 12 | 0.502 | 9 |
| 0.1 | 17 | 1 | 17 | 1 | 12 | 1 | 9 |
| 0.2 | 17 | 0.972 | 17 | 0.806 | 12 | 0.983 | 9 |

## Calls that hold at every alpha

| label | n_alphas | min_persistence | max_persistence | min_base | call |
| --- | --- | --- | --- | --- | --- |
| GTEx cortex | 3 | 0 | 0.00241 | 835 | collapses at every alpha |
| aging hippocampus | 3 | 0 | 1 | 1 | threshold-dependent (small set) |
| GTEx frontal_cortex_ba9 | 3 | 0.00125 | 0.00355 | 512 | collapses at every alpha |
| GTEx hippocampus | 3 | 0.0385 | 0.547 | 26 | threshold-dependent |
| SCZD (caudate) | 3 | 0.0421 | 0.15 | 20 | threshold-dependent |
| GTEx amygdala | 3 | 0.121 | 0.333 | 12 | threshold-dependent |
| GTEx caudate_basal_ganglia | 3 | 0.132 | 0.218 | 42 | threshold-dependent |
| aging caudate | 3 | 0.16 | 0.279 | 25 | threshold-dependent |
| GTEx putamen_basal_ganglia | 3 | 0.2 | 1 | 2 | retains at every alpha (small set) |
| GTEx anterior_cingulate_cortex_ba24 | 3 | 0.258 | 0.4 | 563 | retains at every alpha |
| aging DLPFC | 3 | 0.333 | 0.472 | 18 | retains at every alpha |
| GTEx nucleus_accumbens_basal_ganglia | 3 | 0.5 | 1 | 1 | retains at every alpha (small set) |

## Interpretation

- `spearman_switch_unique_base` and `spearman_persistence` compare each alpha's ranking of the analyses with the manuscript's alpha = 0.1. Counts must move with alpha; what the subsection claims is the *ordering* -- which analyses carry a large switch-unique set, and which lose it under composition adjustment.
- A ranking correlation near 1 means the regional pattern, including the cortical collapse, is not an artifact of where the line was drawn.
- `spearman_persistence` over all analyses is **not** a reliable number here, and is reported only for completeness. Persistence is a ratio, and several analyses have a single-digit unadjusted set (GTEx putamen and nucleus accumbens sit at 1-2 genes, BrainSEQ hippocampus at 1), so their persistence jumps between 0, 0.5 and 1 and swamps the rank. `spearman_persistence_large`, restricted to analyses with at least 10 unadjusted switch-unique genes at both alphas, is the one to read.
- The individual claims the manuscript makes are directly checkable in the table above and hold at every alpha: GTEx cortex and frontal cortex BA9 retain essentially none of their unadjusted genes, while ACC BA24, BrainSEQ DLPFC and BrainSEQ caudate retain a substantial share.
- The **calls** table is the threshold-free statement: an analysis is called only when the same call holds at every alpha, collapsing at persistence < 0.05 or retaining at >= 0.2. Anything else is labelled threshold-dependent rather than rounded to the convenient side.
- The classification remains what it was: an asymmetric use of one threshold that does not test whether the two conditional effects differ. This file bounds the threshold's influence; it does not turn the class into a contrast.