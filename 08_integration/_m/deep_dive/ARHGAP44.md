# ARHGAP44 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Amygdala)
- **Constraint:** LOEUF 0.587  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Amygdala | 0.03 | no | — | no | — | risk allele A decreases ARHGAP44 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Amygdala | rs12936683 | A | -0.2662936747074127 | risk allele A decreases expression of ARHGAP44 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBMS3, RBM41, PPRC1, PABPC4, ENOX1, RBM8A, G3BP2, HNRNPLL, RBM14, AKAP1, MATR3, CELF5, IGF2BP2, G3BP1
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), AGO2(5), AKAP1(5), CELF4(5), CELF5(5), CNOT4(5), CELF6(5), CPEB2(5), DDX19B(5), EIF4A3(5), DHX58(5), G3BP2(5), G3BP1(5), FXR2(5), ESRP2(5), ESRP1(5), ENOX1(5), RBM8A(5), RBMS3(5)

## 5. Interpretation
ARHGAP44 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
