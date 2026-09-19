# NAP1L1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.216  ·  missense o/e 0.60

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | no | yes | — | risk allele T increases usage of junction chr12:76060279-76074203(-) (ENST00000547773.5; n |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs1059143 | T | 0.5517937541007996 | risk allele T increases intron usage of chr12:76060279-76074203(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** DDX58, DHX9, PABPC5, ZNF346, PPRC1, ADAR, IFIH1, RBM41, CNOT4
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), CELF4(3), AKAP1(3), CELF5(3), CELF6(3), CSTF2(3), CNOT4(3), ENOX1(3), ESRP1(3), DDX58(3), DHX58(3), ESRP2(3), GRSF1(3), G3BP1(3), FXR2(3), IGF2BP1(3), IGF2BP2(3), KHDRBS3(3), HNRNPA3(3)

## 5. Interpretation
NAP1L1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
