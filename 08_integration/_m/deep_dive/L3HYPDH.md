# L3HYPDH — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.513  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | no | — | risk allele T increases usage of junction chr14:59473090-59474439(-) (ENST00000463432.5,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele T increases usage of junction chr14:59479351-59505189(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | no | — | risk allele T decreases usage of junction chr14:59504669-59505189(-) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs111541756 | T | 0.9725726246833801 | risk allele T increases intron usage of chr14:59473090-59474439(-) |
| scz | Brain_Cerebellar_Hemisphere | rs111541756 | T | 0.836870014667511 | risk allele T increases intron usage of chr14:59479351-59505189(-) |
| scz | Brain_Cerebellar_Hemisphere | rs111541756 | T | -0.9316191077232361 | risk allele T decreases intron usage of chr14:59504669-59505189(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, RBMS3, RBM41, ELAVL3, AGO1, RBMY1A1, HNRNPLL, NOVA2, PABPC1, MATR3, CELF5, AKAP1, PUM1, PABPC3, DHX9
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), ADAR(4), AGO1(4), AGO2(4), CELF4(4), CELF5(4), CELF6(4), CSTF2(4), CPEB4(4), DAZAP1(4), DDX19B(4), EIF4A3(4), DDX58(4), DHX58(4), DHX9(4), ESRP1(4), ELAVL3(4), GRSF1(4), HNRNPA0(4), NELFE(4)

## 5. Interpretation
L3HYPDH has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
