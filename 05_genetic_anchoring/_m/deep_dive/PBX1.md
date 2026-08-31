# PBX1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cortex)
- **Constraint:** LOEUF 0.129  ·  missense o/e 0.43

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.02 | no | — | yes | — | risk allele T decreases PBX1 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), AGO2(2), CELF4(2), CELF5(2), CELF6(2), FXR1(2), CNOT4(2), DHX58(2), EIF4B(2), ESRP1(2), GRSF1(2), G3BP2(2), FXR2(2), ZFP36L2(2), SNRNP70(2), ZCRB1(2), HNRNPM(2), IGF2BP2(2), IGF2BP1(2), PABPC3(2)

## 5. Interpretation
PBX1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
