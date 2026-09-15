# GRAMD1B — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.683  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele A increases usage of junction chr11:123526171-123540995(+) (unmapped transcrip |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele A decreases usage of junction chr11:123526171-123577367(+) (ENST00000456860.6, |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, G3BP1, SNRPB2, RBMS3, A1CF, AKAP1, RBM14, ZC3H10, ESRP1, YTHDC1, G3BP2, CNOT4, ESRP2, PPRC1, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CNOT4(3), CPEB1(3), CPEB2(3), CSTF2(3), DDX19B(3), EIF4A3(3), DHX58(3), EIF4B(3), ELAVL3(3), ZNF638(3), ESRP1(3), ESRP2(3), FXR1(3)

## 5. Interpretation
GRAMD1B has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
