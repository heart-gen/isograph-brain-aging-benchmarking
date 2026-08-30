# Table 2 -- sQTL/eQTL splicing-specificity contrast

**Table 2. Splicing-QTL are spared relative to expression-QTL specifically in IsoGraph's phenotype-associated, GO-invisible co-switch modules -- an IsoGraph-only method effect.** Paired within-analysis sQTL-odds-ratio / eQTL-odds-ratio contrast (>1 = splicing genetics spared over expression genetics), inverse-variance fixed-effect meta-analysis across brain xQTL analyses (k), with 95% CI, p, and I2 heterogeneity. The contrast removes the shared cis-QTL depletion baseline of constrained network genes. The effect concentrates in the phenotype-associated and GO-invisible module sets, where it is also homogeneous across tissues (I2 = 0.00 for GO-invisible); the GO-visible immune/abundance set is the weakest arm and its nominal significance rests on between-tissue heterogeneity (I2 = 0.68) rather than a consistent effect. The primary internal control is the matched WGCNA baselines, which consume identical switch features and show no effect in any module set. This is the statistical anchor of the genetic-anchoring result (cf. per-gene resolution in Table 3, whose colocalization posteriors are individually modest). Verbatim from 05_genetic_anchoring/_m/qtl_anchoring_meta/qtl_anchoring_meta_contrast.parquet.

| Method | Module set | Analyses (k) | sQTL/eQTL ratio | 95% CI | p (FE) | I2 |
| --- | --- | --- | --- | --- | --- | --- |
| IsoGraph | All modules | 17 | 1.068 | 1.04-1.10 | 1.3e-05 | 0.29 |
| IsoGraph | Phenotype-associated | 12 | 1.163 | 1.10-1.23 | 3.6e-07 | 0.47 |
| IsoGraph | GO-invisible (DTU-without-DGE) | 10 | 1.172 | 1.09-1.26 | 2.3e-05 | 0.0 |
| IsoGraph | GO-visible (immune/abundance control) | 11 | 1.104 | 1.01-1.20 | 0.022 | 0.68 |
| wgcna_switch_only | Phenotype-associated | 11 | 0.977 | 0.90-1.06 | 0.57 | 0.0 |
| wgcna_switch_only | GO-invisible (DTU-without-DGE) | 7 | 1.017 | 0.91-1.13 | 0.76 | 0.36 |
| wgcna_multiplex | Phenotype-associated | 9 | 0.983 | 0.94-1.03 | 0.48 | 0.72 |
| wgcna_multiplex | GO-invisible (DTU-without-DGE) | 7 | 0.993 | 0.91-1.08 | 0.86 | 0.0 |
