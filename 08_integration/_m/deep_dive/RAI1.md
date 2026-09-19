# RAI1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.111  ·  missense o/e 0.82

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | — | no | — | risk allele A increases usage of junction chr17:17793288-17803756(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | no | no | — | risk allele A decreases usage of junction chr17:17798513-17803756(+) (ENST00000353383.6,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs11649804 | A | 0.6095499396324158 | risk allele A increases intron usage of chr17:17793288-17803756(+) |
| scz | Brain_Cerebellar_Hemisphere | rs11649804 | A | -0.6280689835548401 | risk allele A decreases intron usage of chr17:17798513-17803756(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBMS3, PABPC4, RBMY1A1, HNRNPLL, RBM14, RBMS1, MATR3, CELF5, IGF2BP2, G3BP1, FXR1, PABPC3, DHX9, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO2(4), AKAP1(4), CELF1(4), CELF2(4), CELF4(4), CELF5(4), CELF6(4), CMTR1(4), CNOT4(4), CPEB1(4), CPEB4(4), CSTF2(4), DAZAP1(4), DDX19B(4), DHX58(4), DHX9(4), EIF4A3(4)

## 5. Interpretation
RAI1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
