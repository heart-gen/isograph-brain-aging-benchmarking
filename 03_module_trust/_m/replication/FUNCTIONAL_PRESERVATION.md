# Functional preservation of cross-cohort matched modules

Method `isograph`, age model `linear`, 130 matched BrainSEQ<->GTEx module pairs (25 with concordant age effects). Median gene Jaccard 0.0167 — the gene-level overlap these measures are asked to look past.

The null re-pairs each BrainSEQ module with a random GTEx module **from the same gene-count decile**, because all three similarity measures increase with module size. `p_emp` is the size-matched permutation p for the overall matched mean; the Mann-Whitney column asks the separate question of whether the age-concordant pairs are more similar than the other matched pairs.

| measure | matched mean | size-matched null | p_emp | concordant | discordant | diff 95% CI | MWU p |
|---|---|---|---|---|---|---|---|
| `go_jaccard` | 0.0256 | 0.0059 | 0.009901 | 0.1250 | 0.0164 | -0.0303 to 0.3636 | 0.1427 |
| `celltype_r` | 0.2604 | 0.2044 | 0.0396 | 0.3884 | 0.2299 | 0.0528 to 0.2750 | 0.08075 |
| `structure_r` | n/a | n/a | n/a | n/a | n/a | n/a | n/a |

## Verdict

Matched pairs exceed size-matched random pairs on `go_jaccard`, `celltype_r`. The manuscript may state that the cohorts recover related biological programs at a higher level of organisation than gene identity, citing these measures specifically — not as a general claim of replication.
