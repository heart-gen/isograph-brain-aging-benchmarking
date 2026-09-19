# TMEM8B — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Anterior_cingulate_cortex_BA24)
- **Constraint:** LOEUF 0.819  ·  missense o/e 0.86

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr9:35829319-35829877(+) (ENST00000377991.9,ENS |
| als | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.01 | no | no | no | — | risk allele G increases usage of junction chr9:35829391-35829877(+) (ENST00000439587.6; no |
| als | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.01 | no | no | no | — | risk allele G increases usage of junction chr9:35829955-35834461(+) (ENST00000377988.6,ENS |
| als | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr9:35834650-35835011(+) (ENST00000377988.6,ENS |
| als | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.01 | no | — | no | — | risk allele G increases usage of junction chr9:35834650-35841134(+) (unmapped transcript;  |
| als | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr9:35835218-35841134(+) (ENST00000377988.6,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Anterior_cingulate_cortex_BA24 | rs2236291 | G | -0.55522221326828 | risk allele G decreases intron usage of chr9:35829319-35829877(+) |
| als | Brain_Anterior_cingulate_cortex_BA24 | rs2236291 | G | 0.4708060324192047 | risk allele G increases intron usage of chr9:35829391-35829877(+) |
| als | Brain_Anterior_cingulate_cortex_BA24 | rs2236291 | G | 0.7957723736763 | risk allele G increases intron usage of chr9:35829955-35834461(+) |
| als | Brain_Anterior_cingulate_cortex_BA24 | rs2236291 | G | -0.7049988508224487 | risk allele G decreases intron usage of chr9:35834650-35835011(+) |
| als | Brain_Anterior_cingulate_cortex_BA24 | rs2236291 | G | 0.9730847477912903 | risk allele G increases intron usage of chr9:35834650-35841134(+) |
| als | Brain_Anterior_cingulate_cortex_BA24 | rs2236291 | G | -0.6231356263160706 | risk allele G decreases intron usage of chr9:35835218-35841134(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), AGO2(1), ANKHD1(1), CELF6(1), CMTR1(1), CNOT4(1), CPEB1(1), DDX19B(1), DHX58(1), EIF4A3(1), EIF4B(1), ENOX1(1), ESRP1(1), GRSF1(1), HNRNPA0(1), HNRNPA1L2(1), HNRNPA3(1), HNRNPAB(1), HNRNPM(1), HNRNPU(1)

## 5. Interpretation
TMEM8B has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
