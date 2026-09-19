# CACNA1E — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.151  ·  missense o/e 0.46

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele A decreases CACNA1E expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellum | rs7512985 | A | -0.3528773784637451 | risk allele A decreases expression of CACNA1E |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBMS3, RBM41, SNRPB2, RBMS1, PABPC4, ZC3H10, HNRNPLL, RBM14, AKAP1, MATR3, CELF5, IGF2BP2, G3BP1, FXR1
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), AGO2(4), CELF5(4), AKAP1(4), ANKHD1(4), CELF4(4), CNOT4(4), CELF6(4), CPEB2(4), DHX58(4), ZCRB1(4), DDX19B(4), EIF4B(4), EIF4A3(4), ESRP2(4), G3BP1(4), FXR2(4), FXR1(4), HNRNPAB(4), HNRNPCL1(4)

## 5. Interpretation
CACNA1E colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
