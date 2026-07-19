# TTC19 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.049  ·  missense o/e 1.06

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Putamen_basal_ganglia | 0.01 | no | no | yes | — | risk allele T decreases usage of junction chr17:16000245-16001915(+) (ENST00000261647.10,E |
| pd | sQTL | Brain_Putamen_basal_ganglia | 0.01 | no | no | yes | — | risk allele T increases usage of junction chr17:16000382-16001915(+) (ENST00000466729.5,EN |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM42
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), CELF6(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX58(3), DHX58(3), EIF4A3(3), EIF4B(3), ENOX1(3), ESRP1(3), ESRP2(3), FXR2(3), G3BP2(3), HNRNPCL1(3), HNRNPLL(3), IFIH1(3), IGF2BP1(3)

## 5. Interpretation
TTC19 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
