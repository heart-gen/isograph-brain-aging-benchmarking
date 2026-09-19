# STX16 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.961  ·  missense o/e 0.87

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele A decreases usage of junction chr20:58670603-58673631(+) (ENST00000438253.1,EN |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele A increases usage of junction chr20:58673711-58691724(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs218476 | A | -0.49416354298591614 | risk allele A decreases intron usage of chr20:58670603-58673631(+) |
| scz | Brain_Cerebellum | rs218476 | A | 0.34717920422554016 | risk allele A increases intron usage of chr20:58673711-58691724(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM42, A1CF, RBMS3, RBM41, SAMD4A, PPRC1, ADAR, DHX9, HNRNPLL, RBMS1, MATR3, CELF5, IGF2BP2, PABPC5
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ADAR(3), AGO2(3), AKAP1(3), ANKHD1(3), CELF4(3), CELF5(3), CPEB2(3), IFIH1(3), CSTF2(3), DHX58(3), DHX9(3), GRSF1(3), HNRNPA1L2(3), HNRNPLL(3), HNRNPK(3), IGHMBP2(3), IGF2BP2(3), HNRNPM(3), HNRNPU(3)

## 5. Interpretation
STX16 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
