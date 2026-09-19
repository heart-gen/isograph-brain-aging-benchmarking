# FANCL — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible · replicates in BrainSeq

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.348  ·  missense o/e 1.02

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr2:58204226-58221942(-) (ENST00000233741.9,ENS |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | yes | — | risk allele G decreases usage of junction chr2:58222042-58226728(-) (ENST00000233741.9,ENS |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | yes | — | risk allele G decreases usage of junction chr2:58226784-58229814(-) (ENST00000233741.9,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs17049351 | G | 0.8065463304519653 | risk allele G increases intron usage of chr2:58204226-58221942(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs17049351 | G | -1.1809700727462769 | risk allele G decreases intron usage of chr2:58222042-58226728(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs17049351 | G | -1.056998610496521 | risk allele G decreases intron usage of chr2:58226784-58229814(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PABPC4, RBMS3, RBMY1A1, SNRPB2, RBMS1, A1CF, G3BP2, IGF2BP2, FXR1, AKAP1, PABPC5, PPRC1, TARDBP, RBM46, PTBP2
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), CPEB2(5), AGO2(5), ENOX1(5), FXR1(5), DHX58(5), CSTF2(5), RBM28(5), PABPC4(5), MSI1(5), MATR3(5), IGHMBP2(5), KHDRBS2(5), KHDRBS3(5), IFIH1(5), HNRNPLL(5), HNRNPU(5), GRSF1(5), G3BP1(5), SRSF11(5)

## 5. Interpretation
FANCL has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
