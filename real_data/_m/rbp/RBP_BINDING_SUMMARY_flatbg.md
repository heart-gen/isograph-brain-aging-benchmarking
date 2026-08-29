# RBP binding-evidence support for predicted co-switch regulons

- Testable regulon RBPs (sig ∩ ENCODE eCLIP): **130**; peak files present: **39**.

## Per-RBP binding support (independence-respecting headline)

Unit = one **two-sided exact McNemar per RBP** over its **unique nominated regulon genes** (a gene member of a module where the RBP is a significant motif regulon), deduplicated across regions because a gene's switched/constitutive binding is region-invariant; BH-corrected across the 39 testable RBPs. The earlier per-(region×module×RBP) regulon rows are **not** an independent FDR family — the same gene×RBP binding recurs across regions — so they are reported descriptively below, not as the significance count.

- RBPs with preferential switched-interval binding: **31/39**; **binding-supported** (preferential AND BH q≤0.05): **17/39**.
- Median switched−constitutive bound-rate gap across RBPs: **0.013** — small and near-universal. This is binding *capacity* at alternative vs constitutive exons, **not** factor-specific occupancy; it does not establish that the predicted regulon factors selectively bind the switched sequence beyond a generic alternative-exon skew.

### Binding-supported RBPs (switched > constitutive)

| rbp | n_genes | rate_switched | rate_constitutive | rate_diff | mcnemar_fdr |
| --- | --- | --- | --- | --- | --- |
| MATR3 | 2174 | 0.492 | 0.377 | 0.115 | 0.0 |
| HNRNPC | 1036 | 0.642 | 0.571 | 0.07 | 0.0 |
| KHDRBS1 | 3205 | 0.559 | 0.526 | 0.032 | 0.0 |
| QKI | 1830 | 0.548 | 0.506 | 0.042 | 0.0 |
| NONO | 1584 | 0.592 | 0.553 | 0.039 | 0.0 |
| TIAL1 | 1065 | 0.489 | 0.437 | 0.053 | 0.0 |
| HNRNPK | 244 | 0.803 | 0.697 | 0.107 | 0.0 |
| ELAVL1 | 1002 | 0.578 | 0.538 | 0.04 | 0.0003 |
| IGF2BP2 | 2267 | 0.704 | 0.691 | 0.013 | 0.0047 |
| PCBP2 | 263 | 0.285 | 0.232 | 0.053 | 0.0078 |
| FXR1 | 1666 | 0.721 | 0.708 | 0.013 | 0.0078 |
| ZRANB2 | 1091 | 0.545 | 0.519 | 0.027 | 0.0078 |
| IGF2BP3 | 829 | 0.671 | 0.648 | 0.023 | 0.0244 |
| GRSF1 | 406 | 0.618 | 0.569 | 0.049 | 0.0244 |
| SFPQ | 773 | 0.371 | 0.342 | 0.03 | 0.0258 |
| RBFOX2 | 697 | 0.706 | 0.683 | 0.023 | 0.026 |
| AKAP1 | 2057 | 0.656 | 0.641 | 0.015 | 0.0349 |

## Regulon annotation (descriptive, not an FDR family)

- Significant module×RBP regulons with a testable RBP: **192**.
- Regulons whose RBP is binding-supported and that show the preferential direction: **94**; of those GO-invisible: **12**.


_Metric: the RBP's eCLIP peaks are overlapped with each member gene's SWITCHED exon interval (present in one switch isoform) vs its CONSTITUTIVE interval (shared by both); an RBP is scored over its unique nominated regulon genes (two-sided exact McNemar on discordant genes, BH across RBPs). This within-gene null neutralizes peak-dense RBPs. Caveats: ENCODE eCLIP is HepG2/K562, not brain — binding capacity, not brain-specific occupancy; only regulon RBPs with an eCLIP experiment are testable; the switched-vs-constitutive skew is small and shared by most RBPs, so this supports a subset of factors' binding *capacity*, not selective regulon occupancy._