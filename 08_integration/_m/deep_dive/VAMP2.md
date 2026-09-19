# VAMP2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cortex)
- **Constraint:** LOEUF 0.176  ·  missense o/e 0.62

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cortex | 0.01 | yes | yes | no | no annotated structural change | risk allele C decreases usage of junction chr17:8162369-8162878(-) (ENST00000316509.11,ENS |
| als | sQTL | Brain_Cortex | 0.01 | no | — | no | — | risk allele C increases usage of junction chr17:8174272-8174520(-) (unmapped transcript; n |
| als | sQTL | Brain_Cortex | 0.01 | no | — | no | — | risk allele C increases usage of junction chr17:8174272-8176200(-) (unmapped transcript; n |
| als | sQTL | Brain_Cortex | 0.01 | no | — | no | — | risk allele C increases usage of junction chr17:8176026-8176200(-) (unmapped transcript; n |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cortex | rs8066511 | C | -0.38942965865135193 | risk allele C decreases intron usage of chr17:8162369-8162878(-) |
| als | Brain_Cortex | rs8066511 | C | 0.3183293342590332 | risk allele C increases intron usage of chr17:8174272-8174520(-) |
| als | Brain_Cortex | rs8066511 | C | 0.2887972593307495 | risk allele C increases intron usage of chr17:8174272-8176200(-) |
| als | Brain_Cortex | rs8066511 | C | 0.40533390641212463 | risk allele C increases intron usage of chr17:8176026-8176200(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBMS3, PPRC1, PABPC4, RBM14, MATR3, CELF5, IGF2BP2, PABPC3, AKAP1, RBM8A, PTBP2, TARDBP, YBX2
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), AGO1(5), AGO2(5), AKAP1(5), CELF1(5), CELF2(5), CELF4(5), CELF5(5), CELF6(5), CPEB1(5), CPEB2(5), CPEB4(5), CSTF2(5), DAZAP1(5), DDX19B(5), DHX58(5), EIF4A3(5), ELAVL3(5), ESRP1(5), FXR2(5)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for VAMP2 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.
