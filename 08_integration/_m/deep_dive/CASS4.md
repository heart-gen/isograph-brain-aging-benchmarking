# CASS4 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.13 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.523  ·  missense o/e 0.81

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Frontal_Cortex_BA9 | 0.13 | no | — | yes | — | risk allele A decreases CASS4 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Frontal_Cortex_BA9 | rs6014724 | A | -0.4447765350341797 | risk allele A decreases expression of CASS4 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41
- **All switched-motif RBPs (recurrence across regions):** ACO1(1), AGO2(1), AKAP1(1), CSTF2(1), DDX19B(1), DHX58(1), EIF4A3(1), ESRP1(1), FXR1(1), G3BP1(1), GRSF1(1), HNRNPAB(1), HNRNPCL1(1), HNRNPM(1), IGHMBP2(1), KHDRBS2(1), KHDRBS3(1), NELFE(1), PABPC4(1), PABPC5(1)

## 5. Interpretation
CASS4 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
