# TYK2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** LBD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.622  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.03 | no | — | yes | — | risk allele A decreases TYK2 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RALY, HNRNPCL1, A1CF, MATR3, CELF6, IGF2BP1, ENOX1, PABPN1, ESRP2, EIF4A3, CELF5, CELF4, PTBP2
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), EIF4A3(3), CPEB1(3), CELF6(3), SUPV3L1(3), HNRNPA0(3), RBM14(3), ESRP2(3), HNRNPCL1(3), RALY(3), AGO2(2), CELF4(2), CELF5(2), G3BP1(2), HNRNPA3(2), HNRNPAB(2), HNRNPD(2), FXR2(2), CPEB4(2)

## 5. Interpretation
TYK2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
