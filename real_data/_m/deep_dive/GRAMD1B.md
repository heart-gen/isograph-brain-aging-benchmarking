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
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS3, HNRNPCL1, A1CF, ESRP1, G3BP2, ZFP36L2, HNRNPA1L2, RBM24, RALY, SNRNP70, RBM46, DHX58, RBMS1, MSI1, ZNF638
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), CELF6(3), CELF4(3), G3BP2(3), G3BP1(3), CPEB1(3), CSTF2(3), EIF4A3(3), DHX58(3), ENOX1(3), EIF4B(3), HNRNPA0(3), FXR1(3), ESRP1(3), HNRNPLL(3), HNRNPCL1(3), HNRNPAB(3), HNRNPA3(3), HNRNPA1L2(3)

## 5. Interpretation
GRAMD1B has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
