# MEF2C — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.212  ·  missense o/e 0.46

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele G increases MEF2C expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs17558396 | G | 0.28206542134284973 | risk allele G increases expression of MEF2C |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBMS3, SNRPB2, G3BP2, G3BP1, PABPC5
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ESRP2(4), ENOX1(4), DHX58(4), RBFOX2(4), HNRNPA3(4), SNRNP70(4), G3BP2(4), RBM42(4), CPEB2(3), ANKHD1(3), CNOT4(3), PABPC5(3), G3BP1(3), RBM28(3), RBM24(3), RBMS3(3), SNRPB2(3), SRSF4(3), ZFP36L2(3)

## 5. Interpretation
MEF2C colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
