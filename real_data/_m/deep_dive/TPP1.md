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
- **Switched *and* module-enriched (q<0.05) RBPs:** ESRP1, HNRNPA1L2, RBM6
- **All switched-motif RBPs (recurrence across regions):** CELF4(3), ESRP1(3), HNRNPAB(3), HNRNPA1L2(3), ZC3H10(3), RBM6(2), G3BP1(1), CPEB1(1)

## 5. Interpretation
TPP1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
