# PBX1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cortex)
- **Constraint:** LOEUF 0.129  ·  missense o/e 0.43

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.02 | no | — | yes | — | risk allele T decreases PBX1 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, ZCRB1
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), CELF4(1), CELF5(1), CELF6(1), DHX58(1), EIF4B(1), ESRP1(1), G3BP1(1), G3BP2(1), HNRNPA1L2(1), IGF2BP1(1), IGF2BP2(1), MATR3(1), RBFOX2(1), RBM24(1), RBM28(1), RBM41(1), RBM42(1), RBM6(1), RBMS3(1)

## 5. Interpretation
PBX1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
