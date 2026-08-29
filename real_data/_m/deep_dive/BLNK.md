# BLNK — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Amygdala)
- **Constraint:** LOEUF 0.510  ·  missense o/e 0.78

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Amygdala | 0.03 | no | — | no | — | risk allele A increases BLNK expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB1, U2AF2, CPEB4, SYNCRIP, KHDRBS3, RBM4, DDX19B, HNRNPDL, PUM1, ACO1, IFIH1
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), AGO1(5), AGO2(5), AKAP1(5), DDX19B(5), CSTF2(5), ENOX1(5), DHX9(5), G3BP2(5), HNRNPAB(5), ESRP2(5), G3BP1(5), HNRNPDL(5), HNRNPLL(5), IGF2BP1(5), IFIH1(5), ZRANB2(5), SNRPB2(5), KHDRBS2(5), LIN28A(5)

## 5. Interpretation
BLNK colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
