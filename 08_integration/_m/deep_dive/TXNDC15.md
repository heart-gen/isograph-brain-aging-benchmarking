# TXNDC15 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.896  ·  missense o/e 0.79

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Caudate_basal_ganglia | 0.03 | no | no | yes | — | risk allele T decreases usage of junction chr5:134874530-134893492(+) (ENST00000511070.5;  |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ADAR(4), CELF4(4), CELF5(4), CELF6(4), CPEB2(4), CSTF2(4), CPEB4(4), DDX58(4), DDX19B(4), EIF4B(4), ELAVL3(4), DHX9(4), EIF4A3(4), ENOX1(4), GRSF1(4), G3BP1(4), FXR2(4), YTHDC1(4), TRA2A(4)

## 5. Interpretation
TXNDC15 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
