# Table 2 -- sQTL/eQTL splicing-specificity contrast

**Table 2. Splicing-QTL are spared relative to expression-QTL specifically in IsoGraph's phenotype-associated, GO-invisible co-switch modules -- an IsoGraph-only method effect.** Paired within-analysis sQTL-odds-ratio / eQTL-odds-ratio contrast (>1 = splicing genetics spared over expression genetics), inverse-variance fixed-effect meta-analysis across brain xQTL analyses (k), with 95% CI, p, and I2 heterogeneity. The contrast removes the shared cis-QTL depletion baseline of constrained network genes. The effect concentrates in the phenotype-associated and GO-invisible module sets and is null in the GO-visible immune/abundance control; matched WGCNA baselines on identical switch features show no effect. This is the statistical anchor of the genetic-anchoring result (cf. per-gene resolution in Table 3, whose colocalization posteriors are individually modest). Verbatim from real_data/_m/qtl_anchoring_meta/qtl_anchoring_meta_contrast.parquet.

| Method | Module set | Analyses (k) | sQTL/eQTL ratio | 95% CI | p (FE) | I2 |
| --- | --- | --- | --- | --- | --- | --- |
| IsoGraph | All modules | 17 | 1.071 | 1.04-1.10 | 7.8e-06 | 0.35 |
| IsoGraph | Phenotype-associated | 12 | 1.127 | 1.06-1.20 | 0.00015 | 0.22 |
| IsoGraph | GO-invisible (DTU-without-DGE) | 10 | 1.125 | 1.04-1.22 | 0.0026 | 0.15 |
| IsoGraph | GO-visible (immune/abundance control) | 11 | 1.037 | 0.94-1.14 | 0.46 | 0.49 |
| wgcna_switch_only | Phenotype-associated | 11 | 0.977 | 0.90-1.06 | 0.57 | 0.0 |
| wgcna_switch_only | GO-invisible (DTU-without-DGE) | 7 | 1.017 | 0.91-1.13 | 0.76 | 0.36 |
| wgcna_multiplex | Phenotype-associated | 9 | 0.983 | 0.94-1.03 | 0.48 | 0.72 |
| wgcna_multiplex | GO-invisible (DTU-without-DGE) | 7 | 0.993 | 0.91-1.08 | 0.86 | 0.0 |
