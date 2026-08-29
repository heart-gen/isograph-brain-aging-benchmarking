# MED15 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Hypothalamus)
- **Constraint:** LOEUF 0.415  ·  missense o/e 0.80

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Hypothalamus | 0.01 | no | — | yes | — | risk allele T decreases MED15 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), EIF4A3(3), HNRNPCL1(3), IGF2BP2(3), LIN28A(3), PABPC3(3), PABPC5(3), PTBP2(3), PUM2(3), QKI(3), RALY(3), RBM25(3), RBM6(3), SRSF11(3), SRSF4(3), YBX2(3)

## 5. Interpretation
MED15 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
