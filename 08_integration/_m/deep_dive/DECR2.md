# DECR2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF n/a  ·  missense o/e n/a

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele G increases usage of junction chr16:406397-406736(+) (ENST00000631605.1; not i |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr16:406397-407425(+) (ENST00000219481.10,ENST0 |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr16:406777-406834(+) (ENST00000632744.1; not i |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele G increases usage of junction chr16:406833-406890(+) (ENST00000437024.5,ENST00 |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele G increases usage of junction chr16:407223-407425(+) (unmapped transcript; not |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | 0.6874793767929077 | risk allele G increases intron usage of chr16:406397-406736(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | -0.4690103828907013 | risk allele G decreases intron usage of chr16:406397-407425(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | -1.1285455226898193 | risk allele G decreases intron usage of chr16:406777-406834(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | 0.46461474895477295 | risk allele G increases intron usage of chr16:406833-406890(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | 0.8538193702697754 | risk allele G increases intron usage of chr16:407223-407425(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** IGF2BP3, ELAVL3, KHDRBS1, AGO1, RBMY1A1, PABPC1, TRA2A, HNRNPDL
- **All switched-motif RBPs (recurrence across regions):** AGO1(1), AGO2(1), AKAP1(1), CELF4(1), CELF5(1), CMTR1(1), CNOT4(1), CPEB1(1), CPEB4(1), CSTF2(1), DAZAP1(1), DHX58(1), ELAVL3(1), ENOX1(1), G3BP2(1), HNRNPA2B1(1), HNRNPDL(1), HNRNPLL(1), HNRNPM(1), HNRNPU(1)

## 5. Interpretation
DECR2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
