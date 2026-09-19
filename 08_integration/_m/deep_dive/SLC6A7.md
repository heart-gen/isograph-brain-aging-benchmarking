# SLC6A7 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.981  ·  missense o/e 0.90

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele C increases usage of junction chr5:150194911-150196478(+) (unmapped transcript |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele C decreases usage of junction chr5:150194911-150196716(+) (ENST00000230671.7,E |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele C increases usage of junction chr5:150196537-150196716(+) (unmapped transcript |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele C decreases usage of junction chr5:150202450-150202579(+) (ENST00000230671.7,E |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele C increases usage of junction chr5:150202703-150203667(+) (ENST00000230671.7,E |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele C decreases usage of junction chr5:150202703-150203907(+) (unmapped transcript |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele C increases usage of junction chr5:150203779-150203907(+) (ENST00000230671.7,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellar_Hemisphere | rs2240794 | C | 0.6460085511207581 | risk allele C increases intron usage of chr5:150194911-150196478(+) |
| ad | Brain_Cerebellar_Hemisphere | rs2240794 | C | -0.7411422729492188 | risk allele C decreases intron usage of chr5:150194911-150196716(+) |
| ad | Brain_Cerebellar_Hemisphere | rs2240794 | C | 0.5573026537895203 | risk allele C increases intron usage of chr5:150196537-150196716(+) |
| ad | Brain_Cerebellar_Hemisphere | rs2240794 | C | -0.7148227095603943 | risk allele C decreases intron usage of chr5:150202450-150202579(+) |
| ad | Brain_Cerebellar_Hemisphere | rs2240794 | C | 0.6770725846290588 | risk allele C increases intron usage of chr5:150202703-150203667(+) |
| ad | Brain_Cerebellar_Hemisphere | rs2240794 | C | -1.0715936422348022 | risk allele C decreases intron usage of chr5:150202703-150203907(+) |
| ad | Brain_Cerebellar_Hemisphere | rs2240794 | C | 0.46905452013015747 | risk allele C increases intron usage of chr5:150203779-150203907(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SNRPB2, ENOX1, ZC3H10, RBMY1A1, HNRNPLL, RBM14, AKAP1, MATR3, IGF2BP2, G3BP1, FXR1, YTHDC1, TARDBP, HNRNPA3, YBX2
- **All switched-motif RBPs (recurrence across regions):** AGO1(4), AGO2(4), AKAP1(4), CELF6(4), CMTR1(4), CNOT4(4), CPEB1(4), CPEB2(4), CSTF2(4), DHX58(4), EIF4B(4), ELAVL4(4), ENOX1(4), ESRP1(4), ESRP2(4), FXR1(4), FXR2(4), G3BP1(4), HNRNPA1L2(4), HNRNPA3(4)

## 5. Interpretation
SLC6A7 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
