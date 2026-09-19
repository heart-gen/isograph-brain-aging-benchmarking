# Functional preservation of cross-cohort matched modules

Method `isograph`, age model `linear`, 38 matched BrainSEQ<->GTEx module pairs (1 with concordant age effects). Median gene Jaccard 0.0376 — the gene-level overlap these measures are asked to look past.

The null re-pairs each BrainSEQ module with a random GTEx module **from the same gene-count decile**, because all three similarity measures increase with module size. `p_emp` is the size-matched permutation p for the overall matched mean; the Mann-Whitney column asks the separate question of whether the age-concordant pairs are more similar than the other matched pairs.

| measure | n finite | matched mean | size-matched null | p_emp | concordant | discordant | diff 95% CI | MWU p |
|---|---|---|---|---|---|---|---|---|
| `go_jaccard` | 24 | 0.0971 | 0.0296 | 0.000999 | 0.4273 | 0.0828 | n/a | n/a |
| `celltype_r` | 38 | 0.2081 | 0.0396 | 0.005994 | 0.2710 | 0.2064 | n/a | n/a |
| `structure_r` | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a |

## Data coverage

A measure is NaN wherever its upstream analysis was never produced for that region. **These are absent inputs, not null results.**

| measure | regions with no input |
|---|---|
| `structure_r` | brainseq/dlpfc (dlpfc_ba9), gtex/caudate_basal_ganglia (caudate) |

`celltype_composition` is only run for the IsoGraph artifact tree, and `switch_consequence` skips regions with no phenotype-significant switch genes; both are upstream properties, not failures of this analysis. A pair contributes to a measure only when **both** sides have it.

## Verdict

Matched pairs exceed size-matched random pairs on `go_jaccard`, `celltype_r`. The manuscript may state that the cohorts recover related biological programs at a higher level of organisation than gene identity, citing these measures specifically — not as a general claim of replication.

**Untested, not negative:** `structure_r` had no computable pair (see *Data coverage*). The decision rule was applied only to the measures that had data; these remain open.
