# CAMK1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · replicates in BrainSeq

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.973  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | no | — | risk allele G increases usage of junction chr3:9765890-9766125(-) (ENST00000397277.6; not  |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.01 | no | no | no | — | risk allele G increases usage of junction chr3:9766213-9767667(-) (ENST00000397277.6; not  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs62245635 | G | 1.0814597606658936 | risk allele G increases intron usage of chr3:9765890-9766125(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs62245635 | G | 0.9181668758392334 | risk allele G increases intron usage of chr3:9766213-9767667(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBMS3, IGF2BP1, HNRNPA3, CELF4, CELF5, CELF6, MATR3, SNRNP70, A1CF, CPEB2, ZFP36L2, RBM24, AKAP1, ZCRB1
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ADAR(6), AGO1(6), AGO2(6), AKAP1(6), CELF4(6), CELF5(6), CELF6(6), CPEB1(6), CPEB2(6), CPEB4(6), CSTF2(6), DAZAP1(6), DDX19B(6), DHX58(6), DHX9(6), EIF4A3(6), ELAVL1(6), ELAVL3(6), ELAVL4(6)

## 5. Interpretation
CAMK1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
