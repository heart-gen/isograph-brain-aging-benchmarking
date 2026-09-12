# Table 2 -- sQTL/eQTL splicing-specificity contrast

**Table 2. Splicing-QTL are spared relative to expression-QTL in IsoGraph's phenotype-associated co-switch modules -- an IsoGraph-only method effect.** Paired within-analysis sQTL-odds-ratio / eQTL-odds-ratio contrast (>1 = splicing genetics spared over expression genetics), inverse-variance fixed-effect meta-analysis across brain xQTL analyses (k), with 95% CI, p, and I2 heterogeneity. The contrast removes the shared cis-QTL depletion baseline of constrained network genes. The effect is carried by the phenotype-associated set (1.111, p = 3.6e-4); on the 8-tissue set common to all methods it is 1.112 (p = 6.5e-4) at I2 = 0.00. **It does NOT localise to the GO-invisible modules:** GO-invisible (1.068, p = 0.077) and GO-visible (1.084, p = 0.050) are indistinguishable, so this table does not support a GO-invisible-specific genetic claim -- the DTU-without-DGE content claim rests on the GO-invisible gate (Fig S-real-3) instead. The primary internal control is the matched WGCNA baselines, which consume identical switch features and reach significance in no module set; the closest is wgcna_multiplex on the phenotype-associated set (1.046, p = 0.053, I2 = 0.75), so the control is "no baseline clears 0.05", not "every baseline is flat". This is the statistical anchor of the genetic-anchoring result (cf. per-gene resolution in Table 3, whose colocalization posteriors are individually modest). Verbatim from 05_genetic_anchoring/_m/qtl_anchoring_meta/qtl_anchoring_meta_contrast.parquet (regenerated 2026-08-29).

| Method | Module set | Analyses (k) | sQTL/eQTL ratio | 95% CI | p (FE) | I2 |
| --- | --- | --- | --- | --- | --- | --- |
| IsoGraph | All modules | 17 | 1.068 | 1.04-1.10 | 1.3e-05 | 0.29 |
| IsoGraph | Phenotype-associated | 12 | 1.111 | 1.05-1.18 | 0.00036 | 0.23 |
| IsoGraph | GO-invisible (DTU-without-DGE) | 11 | 1.068 | 0.99-1.15 | 0.077 | 0.0 |
| IsoGraph | GO-visible (immune/abundance) | 11 | 1.084 | 1.00-1.17 | 0.05 | 0.28 |
| wgcna_switch_only | Phenotype-associated | 9 | 1.025 | 0.94-1.12 | 0.58 | 0.23 |
| wgcna_switch_only | GO-invisible (DTU-without-DGE) | 8 | 1.066 | 0.98-1.16 | 0.15 | 0.35 |
| wgcna_multiplex | Phenotype-associated | 10 | 1.046 | 1.00-1.10 | 0.053 | 0.75 |
| wgcna_multiplex | GO-invisible (DTU-without-DGE) | 6 | 1.019 | 0.92-1.12 | 0.71 | 0.0 |
