# FAM120AOS — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.176  ·  missense o/e 0.97

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | no | yes | — | risk allele G increases usage of junction chr9:93435922-93437022(-) (ENST00000445280.1; no |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | — | yes | — | risk allele G decreases usage of junction chr9:93437866-93438627(-) (unmapped transcript;  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs7850322 | G | 0.5686156749725342 | risk allele G increases intron usage of chr9:93435922-93437022(-) |
| scz | Brain_Cerebellar_Hemisphere | rs7850322 | G | -0.5794942378997803 | risk allele G decreases intron usage of chr9:93437866-93438627(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ACO1(1), AGO2(1), AKAP1(1), CELF4(1), CELF5(1), CELF6(1), CNOT4(1), CPEB2(1), CSTF2(1), DDX19B(1), EIF4A3(1), ENOX1(1), ESRP1(1), FXR1(1), G3BP2(1), GRSF1(1), HNRNPCL1(1), HNRNPM(1), IGF2BP1(1)

## 5. Interpretation
FAM120AOS has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
