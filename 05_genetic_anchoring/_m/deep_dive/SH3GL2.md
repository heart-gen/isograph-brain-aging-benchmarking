# SH3GL2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.05 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.937  ·  missense o/e 1.01

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Frontal_Cortex_BA9 | 0.05 | no | — | no | — | risk allele T increases SH3GL2 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** TIAL1, ELAVL1, ELAVL4, RC3H1, ZC3H10, CPEB4, RBM14, RBMS1, IGF2BP2, CELF4, SF1
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), AGO1(5), AGO2(5), CELF4(5), CELF5(5), CMTR1(5), CNOT4(5), CPEB1(5), CPEB4(5), CSTF2(5), DAZAP1(5), DDX19B(5), EIF4A3(5), ELAVL1(5), ELAVL2(5), ELAVL3(5), ELAVL4(5), FMR1(5), FXR1(5), GRSF1(5)

## 5. Interpretation
SH3GL2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
