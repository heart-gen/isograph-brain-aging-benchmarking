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
- **Switched *and* module-enriched (q<0.05) RBPs:** KHDRBS1, TIAL1, PPIE, HNRNPC, TIA1, IGF2BP3, HNRNPD, PABPC1, RC3H1, ELAVL3, U2AF2, HNRNPA2B1, DAZAP1, HNRNPA0, HNRNPDL
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), AGO1(5), AGO2(5), AKAP1(5), CPEB2(5), CELF6(5), CPEB1(5), CMTR1(5), CPEB4(5), EIF4B(5), DHX58(5), DDX19B(5), HNRNPC(5), HNRNPAB(5), HNRNPA0(5), EIF4A3(5), ELAVL3(5), ENOX1(5), G3BP1(5)

## 5. Interpretation
PTPRN has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
