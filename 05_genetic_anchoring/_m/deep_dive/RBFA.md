# RBFA — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.977  ·  missense o/e 1.03

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele C increases usage of junction chr18:80042219-80044212(+) (ENST00000306735.10,E |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr18:80048193-80069798(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele C decreases usage of junction chr18:80067308-80069798(+) (unmapped transcript; |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, ZC3H10, RBM14, RBMS3, SNRPB2, A1CF, HNRNPA1L2, PPRC1, HNRNPA3, RBM42, IGF2BP1, RBM28, ENOX1, PABPC3, EIF4B
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO1(3), AGO2(3), AKAP1(3), ANKHD1(3), CELF4(3), CELF5(3), CELF6(3), CMTR1(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), DHX58(3), EIF4A3(3), EIF4B(3), ELAVL3(3), ENOX1(3)

## 5. Interpretation
RBFA has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
