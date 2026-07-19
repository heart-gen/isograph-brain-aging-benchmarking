# BLNK — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Amygdala)
- **Constraint:** LOEUF 0.510  ·  missense o/e 0.78

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Amygdala | 0.03 | no | — | no | — | risk allele A increases BLNK expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB1, CPEB4, KHDRBS3, KHDRBS2, DDX19B, U2AF2, ZRANB2, PUM2, IFIH1, ACO1
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), AGO1(5), AKAP1(5), CSTF2(5), ESRP2(5), ENOX1(5), DHX9(5), G3BP2(5), IGF2BP1(5), IFIH1(5), HNRNPLL(5), PABPC3(5), RBM6(5), RBM25(5), RBM24(5), RBFOX1(5), RBFOX2(5), PABPN1(5), PABPC5(5)

## 5. Interpretation
BLNK colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
