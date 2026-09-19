# Functional preservation of cross-cohort matched modules

Method `wgcna`, age model `linear`, 52 matched BrainSEQ<->GTEx module pairs (2 with concordant age effects). Median gene Jaccard 0.0708 — the gene-level overlap these measures are asked to look past.

The null re-pairs each BrainSEQ module with a random GTEx module **from the same gene-count decile**, because all three similarity measures increase with module size. `p_emp` is the size-matched permutation p for the overall matched mean; the Mann-Whitney column asks the separate question of whether the age-concordant pairs are more similar than the other matched pairs.

| measure | n finite | matched mean | size-matched null | p_emp | concordant | discordant | diff 95% CI | MWU p |
|---|---|---|---|---|---|---|---|---|
| `go_jaccard` | 52 | 0.1356 | 0.1393 | 1 | 0.0857 | 0.1376 | -0.1770 to 0.0711 | n/a |
| `celltype_r` | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| `structure_r` | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a |

## Data coverage

A measure is NaN wherever its upstream analysis was never produced for that region. **These are absent inputs, not null results.**

| measure | regions with no input |
|---|---|
| `celltype_r` | brainseq/caudate (caudate), brainseq/dlpfc (dlpfc_ba9), brainseq/hippocampus (hippocampus), gtex/caudate_basal_ganglia (caudate), gtex/frontal_cortex_ba9 (dlpfc_ba9), gtex/hippocampus (hippocampus) |
| `structure_r` | brainseq/caudate (caudate), brainseq/dlpfc (dlpfc_ba9), brainseq/hippocampus (hippocampus), gtex/caudate_basal_ganglia (caudate), gtex/frontal_cortex_ba9 (dlpfc_ba9), gtex/hippocampus (hippocampus) |

`celltype_composition` is only run for the IsoGraph artifact tree, and `switch_consequence` skips regions with no phenotype-significant switch genes; both are upstream properties, not failures of this analysis. A pair contributes to a measure only when **both** sides have it.

## Verdict

No measure with data exceeds its size-matched null. Per the pre-registered decision rule, the language stays at "matched modules with concordant age effects"; "replicated programs" must not be used.

**Untested, not negative:** `celltype_r`, `structure_r` had no computable pair (see *Data coverage*). The decision rule was applied only to the measures that had data; these remain open.
