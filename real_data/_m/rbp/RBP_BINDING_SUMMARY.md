# RBP binding-evidence support for predicted co-switch regulons

- Testable regulon RBPs (sig ∩ ENCODE eCLIP): **138**; peak files present: **39**.

## Per-RBP binding support (independence-respecting headline)

Unit = one **two-sided exact McNemar per RBP** over its **unique nominated regulon genes** (a gene member of a module where the RBP is a significant motif regulon), deduplicated across regions because a gene's switched/constitutive binding is region-invariant; BH-corrected across the 38 testable RBPs. The earlier per-(region×module×RBP) regulon rows are **not** an independent FDR family — the same gene×RBP binding recurs across regions — so they are reported descriptively below, not as the significance count.

- RBPs with preferential switched-interval binding: **31/38**; **binding-supported** (preferential AND BH q≤0.05): **25/38**, which represent **55** of **66** independent motif families (see `rbp_motif_families.py`).
- Median switched−constitutive bound-rate gap across RBPs: **0.024** — small and near-universal. This is binding *capacity* at alternative vs constitutive exons, **not** factor-specific occupancy; it does not establish that the predicted regulon factors selectively bind the switched sequence beyond a generic alternative-exon skew.

### Binding-supported RBPs (switched > constitutive)

`matched_or` is the Haldane-corrected discordant-pair odds ratio behind the McNemar p; `density_ratio` is peak-covered bp per kb in switched vs constitutive intervals, with a gene-level bootstrap CI. The binary rates alone cannot distinguish a real preference from a longer switched interval having more chance to be touched by some peak.

| rbp | n_genes | rate_switched | rate_constitutive | rate_diff | matched_or | matched_or_ci_low | matched_or_ci_high | density_ratio | density_ci_low | density_ci_high | mcnemar_fdr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MATR3 | 1635 | 0.528 | 0.423 | 0.105 | 6.6393 | 4.5373 | 9.7152 | 1.1487 | 1.0166 | 1.3051 | 0.0 |
| HNRNPC | 2095 | 0.574 | 0.496 | 0.078 | 4.4526 | 3.2505 | 6.0994 | 1.0949 | 1.0161 | 1.1823 | 0.0 |
| TIAL1 | 2067 | 0.565 | 0.505 | 0.06 | 3.5051 | 2.5558 | 4.8068 | 0.8968 | 0.8358 | 0.9647 | 0.0 |
| HNRNPM | 1083 | 0.698 | 0.624 | 0.074 | 4.7209 | 2.9644 | 7.5182 | 1.0133 | 0.9223 | 1.1221 | 0.0 |
| KHDRBS1 | 3188 | 0.576 | 0.539 | 0.037 | 2.5918 | 1.9803 | 3.3923 | 0.7478 | 0.7085 | 0.7877 | 0.0 |
| ELAVL1 | 1771 | 0.645 | 0.591 | 0.054 | 2.8812 | 2.0919 | 3.9682 | 0.8043 | 0.7514 | 0.8622 | 0.0 |
| HNRNPU | 1007 | 0.728 | 0.666 | 0.062 | 4.5429 | 2.7075 | 7.6223 | 0.9425 | 0.8531 | 1.0442 | 0.0 |
| PCBP2 | 406 | 0.655 | 0.539 | 0.116 | 4.2414 | 2.3934 | 7.5163 | 1.0172 | 0.8612 | 1.2108 | 0.0 |
| NONO | 1532 | 0.591 | 0.553 | 0.039 | 2.5325 | 1.7439 | 3.6776 | 0.6639 | 0.6085 | 0.7222 | 0.0 |
| SUPV3L1 | 980 | 0.607 | 0.57 | 0.037 | 3.4828 | 1.9423 | 6.245 | 0.722 | 0.6573 | 0.8008 | 0.0 |
| QKI | 551 | 0.535 | 0.495 | 0.04 | 4.3846 | 1.8704 | 10.2783 | 0.789 | 0.6876 | 0.903 | 0.0007 |
| FXR1 | 1359 | 0.763 | 0.745 | 0.018 | 3.381 | 1.6983 | 6.7307 | 0.6447 | 0.6095 | 0.6792 | 0.0008 |
| HNRNPK | 1067 | 0.76 | 0.724 | 0.037 | 2.0685 | 1.3933 | 3.0708 | 0.818 | 0.7534 | 0.8865 | 0.0008 |
| RBM5 | 1246 | 0.658 | 0.627 | 0.031 | 1.9398 | 1.3338 | 2.8211 | 0.686 | 0.6333 | 0.7425 | 0.0014 |
| PUM2 | 512 | 0.82 | 0.775 | 0.045 | 2.5862 | 1.4107 | 4.7412 | 0.7169 | 0.651 | 0.7858 | 0.0045 |
| IGF2BP2 | 1032 | 0.665 | 0.644 | 0.02 | 2.5556 | 1.3622 | 4.7945 | 0.6258 | 0.5861 | 0.6674 | 0.0073 |
| GRSF1 | 540 | 0.606 | 0.563 | 0.043 | 2.122 | 1.2552 | 3.5873 | 0.6437 | 0.5639 | 0.7361 | 0.0115 |
| TARDBP | 959 | 0.701 | 0.676 | 0.025 | 1.9412 | 1.2039 | 3.1301 | 0.5401 | 0.4993 | 0.5837 | 0.015 |
| IGF2BP3 | 829 | 0.671 | 0.648 | 0.023 | 2.2258 | 1.2224 | 4.0529 | 0.6101 | 0.5625 | 0.6657 | 0.0179 |
| CSTF2 | 503 | 0.843 | 0.811 | 0.032 | 2.6842 | 1.2743 | 5.6541 | 0.849 | 0.7845 | 0.9213 | 0.0179 |
| SFPQ | 498 | 0.468 | 0.428 | 0.04 | 2.0811 | 1.1953 | 3.6232 | 0.7367 | 0.5973 | 0.9182 | 0.0189 |
| MBNL1 | 498 | 0.779 | 0.745 | 0.034 | 2.36 | 1.218 | 4.5728 | 0.5894 | 0.5367 | 0.6504 | 0.0199 |
| PABPC4 | 1810 | 0.666 | 0.651 | 0.015 | 1.5806 | 1.0948 | 2.2821 | 0.6013 | 0.5665 | 0.6394 | 0.0277 |
| RBFOX2 | 657 | 0.7 | 0.679 | 0.021 | 2.3333 | 1.1324 | 4.8078 | 0.6217 | 0.5698 | 0.6791 | 0.0369 |
| CPEB4 | 2092 | 0.405 | 0.393 | 0.012 | 1.4762 | 1.0399 | 2.0955 | 0.513 | 0.4711 | 0.5595 | 0.05 |

## Regulon annotation (descriptive, not an FDR family)

- Significant module×RBP regulons with a testable RBP: **193**.
- Regulons whose RBP is binding-supported and that show the preferential direction: **122**; of those GO-invisible: **24**.


_Metric: the RBP's eCLIP peaks are overlapped with each member gene's SWITCHED exon interval (present in one switch isoform) vs its CONSTITUTIVE interval (shared by both); an RBP is scored over its unique nominated regulon genes (two-sided exact McNemar on discordant genes, BH across RBPs). This within-gene null neutralizes peak-dense RBPs. Caveats: ENCODE eCLIP is HepG2/K562, not brain — binding capacity, not brain-specific occupancy; only regulon RBPs with an eCLIP experiment are testable; the switched-vs-constitutive skew is small and shared by most RBPs, so this supports a subset of factors' binding *capacity*, not selective regulon occupancy._