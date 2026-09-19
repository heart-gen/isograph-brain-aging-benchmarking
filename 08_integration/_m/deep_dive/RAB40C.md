# RAB40C — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.470  ·  missense o/e 0.70

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele A increases usage of junction chr16:590433-625432(+) (unmapped transcript; not |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele A decreases usage of junction chr16:595567-596317(+) (unmapped transcript; not |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele A decreases usage of junction chr16:596448-617208(+) (ENST00000509637.6; not i |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs4984672 | A | 0.3372496962547302 | risk allele A increases intron usage of chr16:590433-625432(+) |
| als | Brain_Cerebellar_Hemisphere | rs4984672 | A | -1.035335659980774 | risk allele A decreases intron usage of chr16:595567-596317(+) |
| als | Brain_Cerebellar_Hemisphere | rs4984672 | A | -0.931286633014679 | risk allele A decreases intron usage of chr16:596448-617208(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, RBMS3, ENOX1, FXR1, MATR3, SNRNP70, YTHDC1, EIF4B, RBMS1, CPEB2, ZFP36L2, IGF2BP2, AKAP1
- **All switched-motif RBPs (recurrence across regions):** AGO1(5), AKAP1(5), DDX19B(5), ELAVL3(5), ENOX1(5), EIF4A3(5), FXR1(5), HNRNPLL(5), HNRNPA1L2(5), G3BP1(5), IGF2BP2(5), NELFE(5), MATR3(5), TRA2A(5), ZRANB2(5), RBMS1(5), RBMS3(5), SUPV3L1(5), SNRNP70(5), SYNCRIP(5)

## 5. Interpretation
RAB40C has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
