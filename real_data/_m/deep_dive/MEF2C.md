# MEF2C — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.212  ·  missense o/e 0.46

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele G increases MEF2C expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** CPEB2(4), G3BP1(4), ENOX1(4), DHX58(4), ZCRB1(4), ZC3H10(4), RBM42(4), PABPC5(4), RBFOX2(4), SNRNP70(3), CELF6(3), HNRNPA1L2(3), CPEB1(1), AKAP1(1), ESRP1(1), MATR3(1), RALY(1), RBM41(1), RBMS1(1), RBMS3(1)

## 5. Interpretation
MEF2C colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
