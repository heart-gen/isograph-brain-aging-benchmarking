# PTPRN — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.735  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellum | 0.02 | no | no | no | — | risk allele G decreases usage of junction chr2:219295141-219296226(-) (ENST00000295718.7,E |
| als | sQTL | Brain_Cerebellum | 0.02 | no | — | no | — | risk allele G increases usage of junction chr2:219295966-219296226(-) (unmapped transcript |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** KHDRBS1, HNRNPD, SYNCRIP, ELAVL3, U2AF2, PABPC1, KHDRBS3, RC3H1, G3BP1, IGHMBP2, KHDRBS2, RBM8A, SRSF10, DDX19B, HNRNPA0
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), AGO1(4), AGO2(4), AKAP1(4), CMTR1(4), CELF6(4), DDX19B(4), CPEB4(4), CPEB1(4), HNRNPCL1(4), G3BP1(4), HNRNPA0(4), DHX58(4), ELAVL3(4), EIF4A3(4), ENOX1(4), ESRP1(4), HNRNPAB(4), HNRNPA3(4)

## 5. Interpretation
PTPRN has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
