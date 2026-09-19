# PTPRU — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.518  ·  missense o/e 0.84

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | no | — | risk allele A decreases PTPRU expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs56335113 | A | -0.19827134907245636 | risk allele A decreases expression of PTPRU |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBMS3, PPRC1, ENOX1, ZC3H10, HNRNPA3, CELF4, CELF5, MATR3, SNRNP70, RBM46, YTHDC1, EIF4B, A1CF, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), AGO1(6), AGO2(6), AKAP1(6), CELF5(6), CELF4(6), CPEB2(6), EIF4A3(6), DDX19B(6), CSTF2(6), ELAVL4(6), ELAVL3(6), EIF4B(6), HNRNPA3(6), HNRNPA0(6), FXR2(6), ESRP2(6), ESRP1(6), ZFP36L2(6)

## 5. Interpretation
PTPRU colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
