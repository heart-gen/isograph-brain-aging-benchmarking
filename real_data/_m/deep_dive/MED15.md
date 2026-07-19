# MED15 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Hypothalamus)
- **Constraint:** LOEUF 0.415  ·  missense o/e 0.80

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Hypothalamus | 0.01 | no | — | yes | — | risk allele T decreases MED15 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZCRB1
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), CELF4(3), CELF5(3), CELF6(3), EIF4A3(3), ENOX1(3), PABPC3(3), PTBP2(3), RALY(3), RBM25(3), RBM6(3), SNRPB2(3), SRSF11(3), YBX2(3), YTHDC1(3), ZC3H10(3), ZCRB1(3)

## 5. Interpretation
MED15 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
