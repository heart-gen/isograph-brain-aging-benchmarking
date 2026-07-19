# TXNDC15 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.896  ·  missense o/e 0.79

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Caudate_basal_ganglia | 0.03 | no | no | yes | — | risk allele T decreases usage of junction chr5:134874530-134893492(+) (ENST00000511070.5;  |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ADAR(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CMTR1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX58(3), DHX9(3), EIF4A3(3), EIF4B(3), ENOX1(3), GRSF1(3), HNRNPAB(3), HNRNPH3(3), HNRNPLL(3), IGF2BP2(3), LIN28A(3)

## 5. Interpretation
TXNDC15 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
