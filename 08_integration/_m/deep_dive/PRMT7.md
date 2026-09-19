# PRMT7 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Hypothalamus)
- **Constraint:** LOEUF 0.935  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Hypothalamus | 0.03 | no | no | yes | — | risk allele T increases usage of junction chr16:68355883-68356701(+) (ENST00000339507.9,EN |
| scz | sQTL | Brain_Hypothalamus | 0.03 | no | no | yes | — | risk allele T decreases usage of junction chr16:68356797-68357054(+) (ENST00000339507.9,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Hypothalamus | rs61733486 | T | 1.3160948753356934 | risk allele T increases intron usage of chr16:68355883-68356701(+) |
| scz | Brain_Hypothalamus | rs61733486 | T | -1.4514796733856201 | risk allele T decreases intron usage of chr16:68356797-68357054(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPC, ELAVL3, KHDRBS1, CPEB1, U2AF2, PPIE
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), AKAP1(3), CELF6(3), CNOT4(3), CPEB1(3), CPEB2(3), DDX19B(3), EIF4A3(3), ELAVL3(3), ESRP1(3), IGHMBP2(3), IGF2BP2(3), HNRNPU(3), HNRNPCL1(3), HNRNPC(3), SUPV3L1(3), SYNCRIP(3), SRSF11(3), U2AF2(3)

## 5. Interpretation
PRMT7 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
