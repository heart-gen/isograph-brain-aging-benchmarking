# THOP1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.844  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele T increases usage of junction chr19:2808444-2810304(+) (ENST00000307741.11,ENS |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr19:2810098-2810304(+) (ENST00000592639.1; not |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr19:2811734-2813115(+) (ENST00000307741.11,ENS |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | no | — | risk allele T increases usage of junction chr19:2812372-2813115(+) (ENST00000395212.8; not |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs941408 | T | 0.48169809579849243 | risk allele T increases intron usage of chr19:2808444-2810304(+) |
| als | Brain_Cerebellar_Hemisphere | rs941408 | T | -0.43667733669281006 | risk allele T decreases intron usage of chr19:2810098-2810304(+) |
| als | Brain_Cerebellar_Hemisphere | rs941408 | T | -0.5076252818107605 | risk allele T decreases intron usage of chr19:2811734-2813115(+) |
| als | Brain_Cerebellar_Hemisphere | rs941408 | T | 0.5851433277130127 | risk allele T increases intron usage of chr19:2812372-2813115(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, RBM14, SNRPB2, ZC3H10, HNRNPLL, CELF5, G3BP1, PABPC3, IGF2BP2, YTHDC1, PTBP2, SNRNP70
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO1(3), AGO2(3), ANKHD1(3), CELF4(3), CELF5(3), CELF6(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), EIF4A3(3), ELAVL1(3), ELAVL3(3), ELAVL4(3), ESRP1(3), FXR2(3), G3BP1(3), GRSF1(3)

## 5. Interpretation
THOP1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
