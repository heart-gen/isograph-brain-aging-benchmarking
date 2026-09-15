# RNASEH2C — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.05 (Brain_Anterior_cingulate_cortex_BA24)
- **Constraint:** LOEUF 1.688  ·  missense o/e 1.28

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Hypothalamus | 0.02 | no | — | yes | — | risk allele A decreases RNASEH2C expression (gene-level; no intron) |
| scz | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.05 | no | no | yes | — | risk allele T increases usage of junction chr11:65715317-65719681(-) (ENST00000644198.1; n |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Hypothalamus | rs7944860 | A | -0.328283429145813 | risk allele A decreases expression of RNASEH2C |
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs4930156 | T | 0.5024814605712891 | risk allele T increases intron usage of chr11:65715317-65719681(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), CELF4(3), CELF5(3), CPEB1(3), CPEB4(3), DHX58(3), SYNCRIP(3), EIF4B(3), ESRP2(3), G3BP1(3), HNRNPCL1(3), HNRNPAB(3), HNRNPU(3), HNRNPM(3), NELFE(3), LIN28A(3), IGHMBP2(3), KHDRBS2(3), PUM1(3)

## 5. Interpretation
RNASEH2C has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
