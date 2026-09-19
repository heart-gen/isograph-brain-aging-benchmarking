# TMEM107 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · replicates in BrainSeq

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cortex)
- **Constraint:** LOEUF 1.472  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cortex | 0.01 | no | — | no | — | risk allele C decreases usage of junction chr17:8162369-8162878(-) (unmapped transcript; n |
| als | sQTL | Brain_Cortex | 0.01 | no | no | no | — | risk allele C increases usage of junction chr17:8174272-8174520(-) (ENST00000316425.9,ENST |
| als | sQTL | Brain_Cortex | 0.01 | no | no | no | — | risk allele C increases usage of junction chr17:8174272-8176200(-) (ENST00000449985.6; not |
| als | sQTL | Brain_Cortex | 0.01 | no | no | no | — | risk allele C increases usage of junction chr17:8176026-8176200(-) (ENST00000316425.9,ENST |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cortex | rs8066511 | C | -0.38942965865135193 | risk allele C decreases intron usage of chr17:8162369-8162878(-) |
| als | Brain_Cortex | rs8066511 | C | 0.3183293342590332 | risk allele C increases intron usage of chr17:8174272-8174520(-) |
| als | Brain_Cortex | rs8066511 | C | 0.2887972593307495 | risk allele C increases intron usage of chr17:8174272-8176200(-) |
| als | Brain_Cortex | rs8066511 | C | 0.40533390641212463 | risk allele C increases intron usage of chr17:8176026-8176200(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS1, PABPC4, ZC3H10, RBM41, AKAP1, RBM6, FXR1
- **All switched-motif RBPs (recurrence across regions):** ACO1(2), AGO1(2), AKAP1(2), CELF6(2), CMTR1(2), CPEB1(2), CPEB2(2), CSTF2(2), DDX19B(2), DHX58(2), EIF4A3(2), EIF4B(2), ESRP1(2), ESRP2(2), FXR1(2), KHDRBS3(2), GRSF1(2), HNRNPCL1(2), HNRNPAB(2), KHDRBS1(2)

## 5. Interpretation
TMEM107 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
