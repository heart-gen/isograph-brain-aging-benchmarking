# SMAP1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cortex)
- **Constraint:** LOEUF 0.452  ·  missense o/e 0.99

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Cortex | 0.02 | no | no | yes | — | risk allele G decreases usage of junction chr6:70858229-70860200(+) (ENST00000316999.9,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Cortex | rs746121 | G | -0.4999880790710449 | risk allele G decreases intron usage of chr6:70858229-70860200(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, RBMS3, ENOX1, CNOT4, IGF2BP1, HNRNPA3, FXR1, CELF6, MATR3, SNRNP70, PABPC5, YTHDC1, RBMS1, SAMD4A, IGF2BP2
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), ADAR(3), AGO1(3), AGO2(3), CELF6(3), CMTR1(3), CNOT4(3), CPEB1(3), CPEB4(3), CSTF2(3), DDX19B(3), DDX58(3), DHX9(3), EIF4A3(3), ELAVL3(3), ENOX1(3), ESRP2(3), FXR1(3), HNRNPA0(3), HNRNPA3(3)

## 5. Interpretation
SMAP1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
