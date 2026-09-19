# FOXN2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.08 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.735  ·  missense o/e 1.10

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Amygdala | 0.02 | no | — | no | — | risk allele G decreases FOXN2 expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellum | 0.08 | no | no | no | — | risk allele G increases usage of junction chr2:48346751-48359047(+) (ENST00000340553.8,ENS |
| scz | sQTL | Brain_Cerebellum | 0.08 | no | no | no | — | risk allele G decreases usage of junction chr2:48359147-48362643(+) (ENST00000340553.8,ENS |
| scz | sQTL | Brain_Cerebellum | 0.08 | no | no | no | — | risk allele G decreases usage of junction chr2:48362707-48373292(+) (ENST00000340553.8,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Amygdala | rs79247094 | G | -0.5428712368011475 | risk allele G decreases expression of FOXN2 |
| scz | Brain_Cerebellum | rs79247094 | G | 1.2604916095733643 | risk allele G increases intron usage of chr2:48346751-48359047(+) |
| scz | Brain_Cerebellum | rs79247094 | G | -0.8725036382675171 | risk allele G decreases intron usage of chr2:48359147-48362643(+) |
| scz | Brain_Cerebellum | rs79247094 | G | -0.8543885350227356 | risk allele G decreases intron usage of chr2:48362707-48373292(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ACO1(1), AGO2(1), AKAP1(1), CPEB1(1), CPEB2(1), CSTF2(1), DDX19B(1), EIF4A3(1), ESRP1(1), G3BP1(1), G3BP2(1), GRSF1(1), HNRNPA1L2(1), HNRNPA3(1), HNRNPAB(1), HNRNPCL1(1), HNRNPLL(1), IGF2BP2(1), KHDRBS2(1)

## 5. Interpretation
FOXN2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
