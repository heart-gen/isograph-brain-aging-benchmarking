# COG7 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.651  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele G decreases COG7 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellum | rs7204605 | G | -0.321371853351593 | risk allele G decreases expression of COG7 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, PABPN1, YBX2, PUM2, CELF6, RBM25, SUPV3L1, KHDRBS1, ZNF638, G3BP2, PABPC1, PPIE, HNRNPLL, FUS, MATR3
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), ADAR(3), AGO2(3), CELF6(3), CMTR1(3), CNOT4(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), DHX9(3), EIF4A3(3), ELAVL3(3), ENOX1(3), FUS(3), G3BP2(3), GRSF1(3), HNRNPA0(3)

## 5. Interpretation
COG7 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
