# FCER1G — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Hypothalamus)
- **Constraint:** LOEUF 0.698  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Hypothalamus | 0.03 | no | — | no | — | risk allele A decreases FCER1G expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB1, KHDRBS3, KHDRBS2, U2AF2, RNASEL, DDX19B, ZRANB2, HNRNPA0, PABPC1, IGHMBP2
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), AKAP1(3), CPEB1(3), DDX19B(3), EIF4A3(3), HNRNPA0(3), HNRNPA3(3), HNRNPCL1(3), HNRNPD(3), HNRNPLL(3), IGF2BP3(3), IGHMBP2(3), KHDRBS2(3), KHDRBS3(3), NELFE(3), NONO(3), NUDT21(3), PABPC1(3), PTBP2(3)

## 5. Interpretation
FCER1G colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
