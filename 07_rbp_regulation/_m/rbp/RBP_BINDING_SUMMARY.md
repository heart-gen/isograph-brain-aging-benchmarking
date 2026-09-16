# RBP binding-evidence support for predicted co-switch regulons

- Testable regulon RBPs (sig ∩ ENCODE eCLIP): **133**; peak files present: **43**.

## Per-RBP binding support (independence-respecting headline)

Unit = one **two-sided exact McNemar per RBP** over its **unique nominated regulon genes** (a gene member of a module where the RBP is a significant motif regulon), deduplicated across regions because a gene's switched/constitutive binding is region-invariant; BH-corrected across the 38 testable RBPs. The earlier per-(region×module×RBP) regulon rows are **not** an independent FDR family — the same gene×RBP binding recurs across regions — so they are reported descriptively below, not as the significance count.

- RBPs with preferential switched-interval binding: **29/38**; **binding-supported** (preferential AND BH q≤0.05): **17/38**, which represent **48** of **73** independent motif families (see `rbp_motif_families.py`).
- Median switched−constitutive bound-rate gap across RBPs: **0.018** — small and near-universal. This is binding *capacity* at alternative vs constitutive exons, **not** factor-specific occupancy; it does not establish that the predicted regulon factors selectively bind the switched sequence beyond a generic alternative-exon skew.

### Binding-supported RBPs (switched > constitutive)

`matched_or` is the Haldane-corrected discordant-pair odds ratio behind the McNemar p; `density_ratio` is peak-covered bp per kb in switched vs constitutive intervals, with a gene-level bootstrap CI. The binary rates alone cannot distinguish a real preference from a longer switched interval having more chance to be touched by some peak.

| rbp | n_genes | rate_switched | rate_constitutive | rate_diff | matched_or | matched_or_ci_low | matched_or_ci_high | density_ratio | density_ci_low | density_ci_high | mcnemar_fdr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| HNRNPC | 4383 | 0.691 | 0.615 | 0.076 | 5.1615 | 4.0656 | 6.5528 | 0.9374 | 0.8984 | 0.9796 | 0.0 |
| MATR3 | 2111 | 0.49 | 0.401 | 0.089 | 5.7595 | 4.1083 | 8.0743 | 1.1487 | 1.0289 | 1.2996 | 0.0 |
| TIAL1 | 3968 | 0.576 | 0.521 | 0.054 | 3.2383 | 2.5775 | 4.0687 | 0.8305 | 0.7874 | 0.8756 | 0.0 |
| ELAVL1 | 3028 | 0.653 | 0.604 | 0.05 | 3.3622 | 2.5407 | 4.4493 | 0.8134 | 0.7742 | 0.8554 | 0.0 |
| KHDRBS1 | 4459 | 0.63 | 0.593 | 0.037 | 2.754 | 2.1737 | 3.4893 | 0.7563 | 0.7254 | 0.7895 | 0.0 |
| HNRNPU | 1365 | 0.834 | 0.765 | 0.069 | 5.3721 | 3.3901 | 8.5129 | 0.8778 | 0.816 | 0.9445 | 0.0 |
| SFPQ | 530 | 0.54 | 0.417 | 0.123 | 9.6667 | 4.5581 | 20.5009 | 1.0149 | 0.8635 | 1.2142 | 0.0 |
| HNRNPK | 1037 | 0.754 | 0.696 | 0.058 | 5.1379 | 2.9273 | 9.018 | 0.885 | 0.8203 | 0.9602 | 0.0 |
| QKI | 1026 | 0.699 | 0.651 | 0.048 | 2.9216 | 1.8635 | 4.5804 | 0.7345 | 0.6509 | 0.8375 | 0.0 |
| PABPC4 | 444 | 0.723 | 0.671 | 0.052 | 6.1111 | 2.2557 | 16.5565 | 0.6535 | 0.5725 | 0.7451 | 0.0001 |
| TIA1 | 1588 | 0.704 | 0.677 | 0.027 | 2.3231 | 1.5399 | 3.5046 | 0.8774 | 0.8188 | 0.9378 | 0.0001 |
| ZRANB2 | 1497 | 0.552 | 0.523 | 0.029 | 2.2394 | 1.5077 | 3.3263 | 0.5929 | 0.5278 | 0.6847 | 0.0001 |
| NONO | 612 | 0.644 | 0.613 | 0.031 | 3.9231 | 1.6581 | 9.2819 | 0.828 | 0.7308 | 0.9278 | 0.0026 |
| PCBP2 | 111 | 0.595 | 0.495 | 0.099 | 23.0 | 1.3554 | 390.3007 | 0.9562 | 0.7134 | 1.2854 | 0.0027 |
| RBM5 | 1106 | 0.699 | 0.675 | 0.024 | 1.8852 | 1.2153 | 2.9244 | 0.6605 | 0.5819 | 0.7781 | 0.0127 |
| RBFOX2 | 2230 | 0.804 | 0.793 | 0.011 | 1.9434 | 1.2164 | 3.105 | 0.6595 | 0.6329 | 0.6846 | 0.0139 |
| GRSF1 | 444 | 0.592 | 0.559 | 0.034 | 3.0 | 1.3129 | 6.8552 | 0.5802 | 0.5004 | 0.675 | 0.0182 |

## Regulon annotation (descriptive, not an FDR family)

- Significant module×RBP regulons with a testable RBP: **123**.
- Regulons whose RBP is binding-supported and that show the preferential direction: **58**; of those GO-invisible: **10**.


_Metric: the RBP's eCLIP peaks are overlapped with each member gene's SWITCHED exon interval (present in one switch isoform) vs its CONSTITUTIVE interval (shared by both); an RBP is scored over its unique nominated regulon genes (two-sided exact McNemar on discordant genes, BH across RBPs). This within-gene null neutralizes peak-dense RBPs. Caveats: ENCODE eCLIP is HepG2/K562, not brain — binding capacity, not brain-specific occupancy; only regulon RBPs with an eCLIP experiment are testable; the switched-vs-constitutive skew is small and shared by most RBPs, so this supports a subset of factors' binding *capacity*, not selective regulon occupancy._