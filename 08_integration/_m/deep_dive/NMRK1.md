# NMRK1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible · replicates in BrainSeq

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.204  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | no | yes | — | risk allele T decreases usage of junction chr9:75069102-75069895(-) (ENST00000376808.8; no |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | no | yes | — | risk allele T increases usage of junction chr9:75069813-75069895(-) (ENST00000361092.9,ENS |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | no | yes | — | risk allele T increases usage of junction chr9:75077580-75078326(-) (ENST00000376811.5; no |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | no | yes | — | risk allele T decreases usage of junction chr9:75083150-75087651(-) (ENST00000466137.1; no |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr9:75069102-75069895(-) (ENST00000376808.8; no |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | no | — | risk allele T increases usage of junction chr9:75069813-75069895(-) (ENST00000361092.9,ENS |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | no | — | risk allele T increases usage of junction chr9:75077580-75078326(-) (ENST00000376811.5; no |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr9:75083150-75087651(-) (ENST00000466137.1; no |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | -0.7377152442932129 | risk allele T decreases intron usage of chr9:75069102-75069895(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | 0.60715651512146 | risk allele T increases intron usage of chr9:75069813-75069895(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | 0.5212721228599548 | risk allele T increases intron usage of chr9:75077580-75078326(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | -0.7121793627738953 | risk allele T decreases intron usage of chr9:75083150-75087651(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | -0.7377152442932129 | risk allele T decreases intron usage of chr9:75069102-75069895(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | 0.60715651512146 | risk allele T increases intron usage of chr9:75069813-75069895(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | 0.5212721228599548 | risk allele T increases intron usage of chr9:75077580-75078326(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs2273766 | T | -0.7121793627738953 | risk allele T decreases intron usage of chr9:75083150-75087651(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, RBM14, RBM42, RBM6, CPEB1, SAMD4A, SRSF10, PABPC1, PPRC1, CPEB4, CNOT4, U2AF2, G3BP2, RBM8A, PPIE
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), AGO1(5), AGO2(5), CELF4(5), CELF5(5), CNOT4(5), CPEB1(5), CPEB4(5), CSTF2(5), DAZAP1(5), DDX58(5), DHX58(5), FXR2(5), G3BP2(5), GRSF1(5), HNRNPA2B1(5), HNRNPA3(5), HNRNPAB(5), HNRNPDL(5), HNRNPK(5)

## 5. Interpretation
NMRK1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
