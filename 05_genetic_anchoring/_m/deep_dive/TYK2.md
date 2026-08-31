# TYK2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** LBD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.622  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.03 | no | — | yes | — | risk allele A decreases TYK2 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBMS3, CELF6, RALY, CNOT4, ZCRB1, MATR3, A1CF, YTHDC1, AKAP1, ESRP2, HNRNPA3, HNRNPCL1, IGF2BP1, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AKAP1(3), GRSF1(3), CELF6(3), CPEB1(3), CNOT4(3), ESRP2(3), HNRNPA0(3), EIF4A3(3), PPIE(3), ZCRB1(3), RBM25(3), RBM14(3), SUPV3L1(3), RBMS3(3), CPEB2(2), CPEB4(2), AGO2(2), HNRNPAB(2)

## 5. Interpretation
TYK2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
