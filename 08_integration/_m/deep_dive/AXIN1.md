# AXIN1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.443  ·  missense o/e 0.97

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele G decreases AXIN1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs11644916 | G | -0.19476677477359772 | risk allele G decreases expression of AXIN1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, PPRC1, IGF2BP3, HNRNPDL, RBM14, RBMS1, CELF5, G3BP1, PABPC3, QKI, AGO1
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CNOT4(3), CPEB1(3), CPEB4(3), CSTF2(3), DAZAP1(3), DDX19B(3), DHX58(3), EIF4B(3), ELAVL1(3), ELAVL2(3), ELAVL3(3), ELAVL4(3)

## 5. Interpretation
AXIN1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
