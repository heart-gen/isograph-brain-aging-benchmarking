# NME4 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.525  ·  missense o/e 1.17

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele G increases usage of junction chr16:406397-406736(+) (unmapped transcript; not |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele G decreases usage of junction chr16:406397-407425(+) (unmapped transcript; not |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele G decreases usage of junction chr16:406777-406834(+) (unmapped transcript; not |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele G increases usage of junction chr16:406833-406890(+) (unmapped transcript; not |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele G increases usage of junction chr16:407223-407425(+) (unmapped transcript; not |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | 0.6874793767929077 | risk allele G increases intron usage of chr16:406397-406736(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | -0.4690103828907013 | risk allele G decreases intron usage of chr16:406397-407425(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | -1.1285455226898193 | risk allele G decreases intron usage of chr16:406777-406834(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | 0.46461474895477295 | risk allele G increases intron usage of chr16:406833-406890(+) |
| als | Brain_Cerebellar_Hemisphere | rs3743890 | G | 0.8538193702697754 | risk allele G increases intron usage of chr16:407223-407425(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, PABPC5
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), AGO1(2), CELF4(2), CELF6(2), CELF5(2), CPEB1(2), CNOT4(2), CSTF2(2), DDX19B(2), CPEB2(2), CPEB4(2), EIF4A3(2), ELAVL3(2), ESRP1(2), ELAVL4(2), ZFP36(2), TARDBP(2), GRSF1(2), HNRNPA0(2), HNRNPCL1(2)

## 5. Interpretation
NME4 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
