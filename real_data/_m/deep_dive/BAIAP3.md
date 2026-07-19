# BAIAP3 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Hypothalamus)
- **Constraint:** LOEUF 1.383  ·  missense o/e 1.17

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Hypothalamus | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr16:1333749-1338540(+) (ENST00000397488.6,ENST |
| als | sQTL | Brain_Hypothalamus | 0.01 | no | no | yes | — | risk allele G decreases usage of junction chr16:1334756-1338540(+) (ENST00000324385.9,ENST |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM8A
- **All switched-motif RBPs (recurrence across regions):** CPEB4(3), ENOX1(3), G3BP1(3), RBM14(3), RBM8A(3), SNRPA(3), YTHDC1(3), ZC3H10(3)

## 5. Interpretation
BAIAP3 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
