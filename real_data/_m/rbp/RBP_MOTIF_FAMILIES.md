# RBP motif families

1194 human ATtRACT matrices for 160 RBPs collapse to **136 similarity families** at a distance cut of 0.25 (mean column correlation >= 0.75 under average linkage, minimum 4 overlapping columns, sense orientation only). 107 families span more than one RBP; the largest holds 112 matrices.

Report motif recurrence **per family**. Counting matrices, or even distinct RBPs, treats near-duplicate motifs as independent evidence.

## Cut calibration

Similarity quantiles for pairs of matrices annotated to the same group vs different groups. The cut is only meaningful if these separate.

| grouping | pairs within | pairs between | within q25 | within q50 | between q75 | between q90 |
|---|---|---|---|---|---|---|
| same RBP | 16,010 | 696,211 | 0.382 | 0.661 | 0.467 | 0.667 |
| same ATtRACT family | 220,560 | 491,661 | 0.111 | 0.333 | 0.467 | 0.667 |

The two distributions overlap heavily: matrices for *different* RBPs are about as similar as matrices for the same RBP. That is the redundancy this analysis exists to absorb — many RBPs bind near-identical short degenerate elements — but it also means no cut is uniquely correct, so the sweep below is reported rather than a single number defended.

## Cut sensitivity

| cut | families | largest | singletons | RBP purity |
|---|---|---|---|---|
| 0.10 | 262 | 92 | 78 | 0.414 |
| 0.15 | 213 | 98 | 48 | 0.385 |
| 0.20 | 172 | 112 | 32 | 0.360 |
| 0.25 **(pinned)** | 136 | 112 | 21 | 0.340 |
| 0.30 | 100 | 140 | 9 | 0.319 |
| 0.35 | 72 | 185 | 4 | 0.297 |
| 0.40 | 57 | 202 | 4 | 0.286 |

## Largest families

| family | label | matrices | RBPs |
|---|---|---|---|
| F0006 | HNRNPH2 | 112 | 20 |
| F0114 | ELAVL1 | 88 | 36 |
| F0025 | SRSF1 | 54 | 19 |
| F0112 | ELAVL4 | 47 | 19 |
| F0028 | SRSF1 | 45 | 21 |
| F0007 | HNRNPA1 | 42 | 13 |
| F0054 | PABPC1 | 38 | 27 |
| F0128 | CELF1 | 35 | 14 |
| F0086 | HNRNPK | 28 | 10 |
| F0105 | PTBP1 | 24 | 6 |
| F0106 | PTBP1 | 24 | 5 |
| F0087 | PCBP2 | 20 | 6 |
| F0067 | HNRNPA1 | 19 | 7 |
| F0080 | PUM1 | 18 | 7 |
| F0107 | PTBP1 | 17 | 6 |
