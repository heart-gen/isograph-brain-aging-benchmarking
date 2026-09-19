# RBP binding-evidence support for predicted co-switch regulons

- Testable regulon RBPs (sig ∩ ENCODE eCLIP): **125**; peak files present: **43**.

## Per-RBP binding support (independence-respecting headline)

Unit = one **two-sided exact McNemar per RBP** over its **unique nominated regulon genes** (a gene member of a module where the RBP is a significant motif regulon), deduplicated across regions because a gene's switched/constitutive binding is region-invariant; BH-corrected across the 35 testable RBPs. The earlier per-(region×module×RBP) regulon rows are **not** an independent FDR family — the same gene×RBP binding recurs across regions — so they are reported descriptively below, not as the significance count.

- RBPs with preferential switched-interval binding: **31/35**; **binding-supported** (preferential AND BH q≤0.05): **21/35**, which represent **47** of **60** independent motif families (see `rbp_motif_families.py`).
- Median switched−constitutive bound-rate gap across RBPs: **0.020** — small and near-universal. This is binding *capacity* at alternative vs constitutive exons, **not** factor-specific occupancy; it does not establish that the predicted regulon factors selectively bind the switched sequence beyond a generic alternative-exon skew.

### Binding-supported RBPs (switched > constitutive)

`matched_or` is the Haldane-corrected discordant-pair odds ratio behind the McNemar p; `density_ratio` is peak-covered bp per kb in switched vs constitutive intervals, with a gene-level bootstrap CI. The binary rates alone cannot distinguish a real preference from a longer switched interval having more chance to be touched by some peak.

| rbp | n_genes | rate_switched | rate_constitutive | rate_diff | matched_or | matched_or_ci_low | matched_or_ci_high | density_ratio | density_ci_low | density_ci_high | mcnemar_fdr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MATR3 | 5001 | 0.448 | 0.342 | 0.106 | 6.8142 | 5.4717 | 8.4861 | 1.1571 | 1.0701 | 1.2585 | 0.0 |
| HNRNPC | 1919 | 0.678 | 0.568 | 0.11 | 6.9718 | 4.9044 | 9.9109 | 1.0048 | 0.9302 | 1.0869 | 0.0 |
| HNRNPM | 2700 | 0.657 | 0.599 | 0.058 | 4.9494 | 3.5161 | 6.9669 | 0.9525 | 0.8947 | 1.0124 | 0.0 |
| TIAL1 | 3087 | 0.592 | 0.537 | 0.055 | 3.1384 | 2.4382 | 4.0395 | 0.8386 | 0.7797 | 0.9008 | 0.0 |
| ELAVL1 | 3338 | 0.7 | 0.66 | 0.04 | 2.6196 | 2.0296 | 3.3812 | 0.8086 | 0.7727 | 0.8459 | 0.0 |
| KHDRBS1 | 3757 | 0.633 | 0.599 | 0.035 | 2.6273 | 2.0325 | 3.3962 | 0.765 | 0.7285 | 0.8054 | 0.0 |
| HNRNPU | 1026 | 0.745 | 0.666 | 0.079 | 4.9512 | 3.0804 | 7.9584 | 0.9397 | 0.842 | 1.058 | 0.0 |
| HNRNPK | 1087 | 0.67 | 0.6 | 0.07 | 5.1081 | 3.1035 | 8.4075 | 0.8194 | 0.7357 | 0.907 | 0.0 |
| QKI | 1392 | 0.624 | 0.563 | 0.061 | 3.4638 | 2.3716 | 5.059 | 0.8253 | 0.7495 | 0.9127 | 0.0 |
| FUS | 3330 | 0.746 | 0.724 | 0.021 | 2.5604 | 1.8176 | 3.6068 | 0.71 | 0.6854 | 0.7359 | 0.0 |
| TARDBP | 2384 | 0.635 | 0.605 | 0.03 | 2.152 | 1.5942 | 2.905 | 0.6323 | 0.5873 | 0.6841 | 0.0 |
| ZRANB2 | 3062 | 0.573 | 0.549 | 0.024 | 1.9172 | 1.4593 | 2.5187 | 0.655 | 0.6119 | 0.7043 | 0.0 |
| NONO | 1587 | 0.618 | 0.588 | 0.03 | 2.6271 | 1.7192 | 4.0144 | 0.8786 | 0.8122 | 0.9552 | 0.0 |
| IGF2BP3 | 2755 | 0.713 | 0.693 | 0.02 | 2.2043 | 1.5587 | 3.1172 | 0.5892 | 0.5541 | 0.6347 | 0.0 |
| AKAP1 | 4406 | 0.691 | 0.673 | 0.017 | 1.7343 | 1.3617 | 2.2089 | 0.663 | 0.6377 | 0.6889 | 0.0 |
| RBM5 | 1540 | 0.64 | 0.608 | 0.032 | 2.2099 | 1.5246 | 3.2031 | 0.6972 | 0.6237 | 0.794 | 0.0 |
| U2AF2 | 2066 | 0.778 | 0.758 | 0.02 | 2.4909 | 1.6003 | 3.8772 | 0.6141 | 0.5839 | 0.6505 | 0.0001 |
| PUM2 | 3556 | 0.659 | 0.641 | 0.017 | 1.7168 | 1.317 | 2.2379 | 0.6612 | 0.6331 | 0.6898 | 0.0001 |
| IGF2BP2 | 5117 | 0.68 | 0.671 | 0.009 | 1.5732 | 1.1856 | 2.0877 | 0.5778 | 0.5607 | 0.5961 | 0.0034 |
| PABPC4 | 5234 | 0.639 | 0.63 | 0.009 | 1.3254 | 1.0703 | 1.6413 | 0.5298 | 0.511 | 0.55 | 0.0191 |
| PCBP1 | 604 | 0.637 | 0.608 | 0.03 | 2.1613 | 1.1837 | 3.9463 | 0.4124 | 0.3599 | 0.4694 | 0.0221 |

## Regulon annotation (descriptive, not an FDR family)

- Significant module×RBP regulons with a testable RBP: **111**.
- Regulons whose RBP is binding-supported and that show the preferential direction: **69**; of those GO-invisible: **10**.


_Metric: the RBP's eCLIP peaks are overlapped with each member gene's SWITCHED exon interval (present in one switch isoform) vs its CONSTITUTIVE interval (shared by both); an RBP is scored over its unique nominated regulon genes (two-sided exact McNemar on discordant genes, BH across RBPs). This within-gene null neutralizes peak-dense RBPs. Caveats: ENCODE eCLIP is HepG2/K562, not brain — binding capacity, not brain-specific occupancy; only regulon RBPs with an eCLIP experiment are testable; the switched-vs-constitutive skew is small and shared by most RBPs, so this supports a subset of factors' binding *capacity*, not selective regulon occupancy._