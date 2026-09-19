# INTS8 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.446  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.03 | no | no | yes | — | risk allele G increases usage of junction chr8:94873477-94874552(+) (ENST00000343161.8,ENS |
| ad | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.03 | no | no | yes | — | risk allele G decreases usage of junction chr8:94874602-94876074(+) (ENST00000343161.8,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Nucleus_accumbens_basal_ganglia | rs13248607 | G | 0.5397234559059143 | risk allele G increases intron usage of chr8:94873477-94874552(+) |
| ad | Brain_Nucleus_accumbens_basal_ganglia | rs13248607 | G | -0.4498647451400757 | risk allele G decreases intron usage of chr8:94874602-94876074(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PUM1, NUDT21, SF1, CMTR1, RBM4
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO2(4), CMTR1(4), CNOT4(4), CPEB2(4), DDX19B(4), DDX58(4), DHX58(4), DHX9(4), EIF4A3(4), EIF4B(4), ENOX1(4), ERI1(4), ESRP1(4), FMR1(4), FXR2(4), G3BP2(4), GRSF1(4)

## 5. Interpretation
INTS8 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
