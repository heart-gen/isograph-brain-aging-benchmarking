# SH3GL2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.05 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.937  ·  missense o/e 1.01

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Frontal_Cortex_BA9 | 0.05 | no | — | no | — | risk allele T increases SH3GL2 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ELAVL3, KHDRBS3, RC3H1, KHDRBS1, HNRNPA2B1, KHDRBS2, DDX19B, U2AF2, AKAP1, AGO2, CPEB4, CELF4, HNRNPA0, SFPQ, PUM2
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CMTR1(3), CPEB1(3), CPEB4(3), CSTF2(3), DDX19B(3), EIF4A3(3), ELAVL2(3), ELAVL3(3), FXR1(3), HNRNPA0(3), HNRNPA1L2(3), HNRNPA2B1(3), HNRNPAB(3)

## 5. Interpretation
SH3GL2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
