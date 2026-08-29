# TPP1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.48 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.746  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellum | 0.48 | no | — | yes | — | risk allele G decreases usage of junction chr11:6611818-6611958(-) (unmapped transcript; n |
| als | sQTL | Brain_Cerebellum | 0.48 | no | no | yes | — | risk allele G increases usage of junction chr11:6614686-6614866(-) (ENST00000299427.12,ENS |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, SNRPB2, RBMS3, AKAP1, RBM6, ZC3H10, ESRP1, YTHDC1, G3BP2, PABPC5, SART3, HNRNPLL, ZCRB1, KHDRBS2
- **All switched-motif RBPs (recurrence across regions):** AKAP1(4), CELF6(4), HNRNPAB(4), ESRP1(4), G3BP2(4), FXR2(4), NELFE(4), LIN28A(4), KHDRBS2(4), HNRNPLL(4), TARDBP(4), YTHDC1(4), SNRPB2(4), SRSF4(4), RBM6(4), RBM25(4), SART3(4), RBMS3(4), ZC3H10(4), ZCRB1(4)

## 5. Interpretation
TPP1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
