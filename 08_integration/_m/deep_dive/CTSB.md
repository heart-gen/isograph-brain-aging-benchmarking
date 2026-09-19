# CTSB — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Amygdala)
- **Constraint:** LOEUF 1.705  ·  missense o/e 1.63

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Amygdala | 0.01 | no | no | no | — | risk allele C increases usage of junction chr8:11844241-11845082(-) (ENST00000530640.7,ENS |
| pd | sQTL | Brain_Amygdala | 0.01 | yes | yes | no | no annotated structural change | risk allele C decreases usage of junction chr8:11845222-11845661(-) (ENST00000345125.8,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Amygdala | rs1293295 | C | 1.077845573425293 | risk allele C increases intron usage of chr8:11844241-11845082(-) |
| pd | Brain_Amygdala | rs1293295 | C | -0.6373957395553589 | risk allele C decreases intron usage of chr8:11845222-11845661(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM6, IGF2BP1, ZC3H10, FXR1, YTHDC1, A1CF, SAMD4A, IGF2BP2, PABPC4, RALY, SART3, PABPN1, NELFE, DHX58, LIN28A
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), AGO2(5), DHX58(5), ELAVL3(5), HNRNPA0(5), HNRNPAB(5), KHDRBS3(5), IGF2BP1(5), KHDRBS2(5), IGHMBP2(5), RBM6(5), RALY(5), PABPC4(5), NELFE(5), ZC3H10(5), SYNCRIP(5), YTHDC1(5), TRA2A(5), SART3(5), SAMD4A(5)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for CTSB in PD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.
