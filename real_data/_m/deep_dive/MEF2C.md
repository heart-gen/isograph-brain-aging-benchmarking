# MEF2C — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.212  ·  missense o/e 0.46

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele G increases MEF2C expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, RBMS1, FXR2, RBMS3, G3BP1, SNRPB2, HNRNPA3, RBM42, RBM28, ENOX1, PABPN1, CELF6, AKAP1, MATR3
- **All switched-motif RBPs (recurrence across regions):** ENOX1(6), HNRNPA3(6), DHX58(6), G3BP2(6), ESRP2(6), RBFOX2(6), SUPV3L1(6), RBM28(6), ZFP36L2(5), G3BP1(5), RBM42(5), SNRPB2(5), RBMS3(5), PABPC5(5), CPEB2(5), ANKHD1(5), CNOT4(5), SRSF4(5), SNRNP70(4), PABPN1(4)

## 5. Interpretation
MEF2C colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
