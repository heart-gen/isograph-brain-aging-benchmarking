# Table 2 -- sQTL/eQTL splicing-specificity contrast

**Table 2. Splicing-QTL are spared relative to expression-QTL in IsoGraph's phenotype-associated co-switch modules -- an IsoGraph-only method effect.** Paired within-analysis sQTL-odds-ratio / eQTL-odds-ratio contrast (>1 = splicing genetics spared over expression genetics), inverse-variance fixed-effect meta-analysis across brain xQTL analyses (k), with 95% CI, p, and I2 heterogeneity. The contrast removes the shared cis-QTL depletion baseline of constrained network genes. The effect is carried by the phenotype-associated set (1.111, p = 3.6e-4); on the 8-tissue set common to all methods it is 1.108 (p = 0.001) at I2 = 0.00. **It does NOT localise to the GO-invisible modules:** GO-invisible (1.068, p = 0.077) and GO-visible (1.084, p = 0.050) are indistinguishable, so this table does not support a GO-invisible-specific genetic claim -- the DTU-without-DGE content claim rests on the GO-invisible gate (Fig S-real-3) instead. The primary internal control is the matched WGCNA baselines, which consume identical switch features and are null in every module set (p >= 0.41). This is the statistical anchor of the genetic-anchoring result (cf. per-gene resolution in Table 3, whose colocalization posteriors are individually modest). Verbatim from 05_genetic_anchoring/_m/qtl_anchoring_meta/qtl_anchoring_meta_contrast.parquet (regenerated 2026-08-29).

| Method | Module set | Analyses (k) | sQTL/eQTL ratio | 95% CI | p (FE) | I2 |
| --- | --- | --- | --- | --- | --- | --- |
| IsoGraph | All modules | 17 | 1.068 | 1.04-1.10 | 1.3e-05 | 0.29 |
| IsoGraph | Phenotype-associated | 12 | 1.111 | 1.05-1.18 | 0.00036 | 0.23 |
| IsoGraph | GO-invisible (DTU-without-DGE) | 11 | 1.068 | 0.99-1.15 | 0.077 | 0.0 |
| IsoGraph | GO-visible (immune/abundance) | 11 | 1.084 | 1.00-1.17 | 0.05 | 0.28 |
| wgcna_switch_only | Phenotype-associated | 11 | 0.977 | 0.90-1.06 | 0.57 | 0.0 |
| wgcna_switch_only | GO-invisible (DTU-without-DGE) | 7 | 1.017 | 0.91-1.13 | 0.76 | 0.36 |
| wgcna_multiplex | Phenotype-associated | 9 | 0.983 | 0.94-1.03 | 0.48 | 0.72 |
| wgcna_multiplex | GO-invisible (DTU-without-DGE) | 7 | 0.993 | 0.91-1.08 | 0.86 | 0.0 |
