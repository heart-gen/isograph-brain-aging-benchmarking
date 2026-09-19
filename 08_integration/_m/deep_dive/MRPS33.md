# MRPS33 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.550  ·  missense o/e 0.80

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.01 | no | no | yes | — | risk allele A increases usage of junction chr7:141006535-141010419(-) (ENST00000324787.10, |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.01 | no | no | yes | — | risk allele A decreases usage of junction chr7:141006535-141014551(-) (ENST00000484502.1;  |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.01 | no | — | yes | — | risk allele A decreases usage of junction chr7:141006535-141014911(-) (unmapped transcript |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Putamen_basal_ganglia | rs151250101 | A | 1.4888312816619873 | risk allele A increases intron usage of chr7:141006535-141010419(-) |
| scz | Brain_Putamen_basal_ganglia | rs151250101 | A | -1.356309175491333 | risk allele A decreases intron usage of chr7:141006535-141014551(-) |
| scz | Brain_Putamen_basal_ganglia | rs151250101 | A | -1.285607099533081 | risk allele A decreases intron usage of chr7:141006535-141014911(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** DHX9, PPRC1, ADAR, RBM41, ERI1, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ACO1(1), ADAR(1), AGO1(1), AGO2(1), CELF4(1), CELF5(1), CELF6(1), CMTR1(1), CPEB2(1), CPEB4(1), CSTF2(1), DAZAP1(1), DDX19B(1), DHX58(1), DHX9(1), EIF4A3(1), EIF4B(1), ERI1(1), ESRP1(1)

## 5. Interpretation
MRPS33 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
