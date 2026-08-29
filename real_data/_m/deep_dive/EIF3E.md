# EIF3E — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.756  ·  missense o/e 0.61

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A decreases usage of junction chr8:108229195-108234998(-) (ENST00000220849.10, |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A increases usage of junction chr8:108233672-108234998(-) (ENST00000521440.6,E |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A increases usage of junction chr8:108241913-108243293(-) (ENST00000518345.2,E |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A decreases usage of junction chr8:108241913-108248613(-) (ENST00000220849.10, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs13257548 | A | -0.4141279458999634 | risk allele A decreases intron usage of chr8:108229195-108234998(-) |
| scz | Brain_Cerebellar_Hemisphere | rs13257548 | A | 0.38281089067459106 | risk allele A increases intron usage of chr8:108233672-108234998(-) |
| scz | Brain_Cerebellar_Hemisphere | rs13257548 | A | 0.3942139446735382 | risk allele A increases intron usage of chr8:108241913-108243293(-) |
| scz | Brain_Cerebellar_Hemisphere | rs13257548 | A | -0.41363778710365295 | risk allele A decreases intron usage of chr8:108241913-108248613(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ADAR(4), AGO2(4), CPEB1(4), DDX19B(4), EIF4A3(4), ELAVL3(4), ESRP1(4), HNRNPA0(4), HNRNPCL1(4), HNRNPM(4), IFIH1(4), IGHMBP2(4), KHDRBS2(4), LIN28A(4), NELFE(4), PABPC4(4), PHAX(4), PTBP2(4), RBFOX1(4), RBM25(4)

## 5. Interpretation
EIF3E has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
