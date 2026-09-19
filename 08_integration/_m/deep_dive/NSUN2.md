# NSUN2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.965  ·  missense o/e 0.99

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele G increases usage of junction chr5:6600232-6602374(-) (unmapped transcript; no |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | yes | — | risk allele G decreases usage of junction chr5:6600232-6602461(-) (ENST00000264670.11,ENST |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr5:6602500-6604138(-) (ENST00000264670.11,ENST |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr5:6605408-6606820(-) (ENST00000264670.11,ENST |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele G increases usage of junction chr5:6600232-6602374(-) (unmapped transcript; no |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr5:6600232-6602461(-) (ENST00000264670.11,ENST |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele G increases usage of junction chr5:6602500-6604138(-) (ENST00000264670.11,ENST |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele G increases usage of junction chr5:6605408-6606820(-) (ENST00000264670.11,ENST |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | 0.4119894504547119 | risk allele G increases intron usage of chr5:6600232-6602374(-) |
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | -0.47249501943588257 | risk allele G decreases intron usage of chr5:6600232-6602461(-) |
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | 0.4680074155330658 | risk allele G increases intron usage of chr5:6602500-6604138(-) |
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | 0.6041921377182007 | risk allele G increases intron usage of chr5:6605408-6606820(-) |
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | 0.4119894504547119 | risk allele G increases intron usage of chr5:6600232-6602374(-) |
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | -0.47249501943588257 | risk allele G decreases intron usage of chr5:6600232-6602461(-) |
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | 0.4680074155330658 | risk allele G increases intron usage of chr5:6602500-6604138(-) |
| scz | Brain_Cerebellar_Hemisphere | rs3756429 | G | 0.6041921377182007 | risk allele G increases intron usage of chr5:6605408-6606820(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, ZFP36L2, A1CF, RBM6, RBMS3, RBM41, SAMD4A, PUM2, RBM8A, ADAR, DHX9, PABPC1, PPIE, HNRNPLL, DDX58
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ADAR(4), AGO1(4), AGO2(4), CPEB2(4), CSTF2(4), DDX19B(4), DDX58(4), DHX58(4), DHX9(4), EIF4A3(4), ENOX1(4), ESRP1(4), FXR1(4), HNRNPA3(4), HNRNPAB(4), HNRNPCL1(4), HNRNPDL(4), HNRNPLL(4), HNRNPU(4)

## 5. Interpretation
NSUN2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
