# CHRNB2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cortex)
- **Constraint:** LOEUF 1.116  ·  missense o/e 0.89

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.01 | no | — | no | — | risk allele C decreases CHRNB2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs11264227 | C | -0.12016735970973969 | risk allele C decreases expression of CHRNB2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM6, G3BP2, ENOX1, IGF2BP1, HNRNPA3, MATR3, PABPC5, YTHDC1, SAMD4A, RBM24, DHX9, FXR2, RALY, SART3
- **All switched-motif RBPs (recurrence across regions):** ADAR(6), AGO1(6), AGO2(6), CELF1(6), CELF2(6), CPEB1(6), CPEB4(6), CSTF2(6), DAZAP1(6), DDX19B(6), DDX58(6), DHX58(6), DHX9(6), EIF4A3(6), ELAVL1(6), ELAVL2(6), ELAVL3(6), ELAVL4(6), ENOX1(6), ESRP1(6)

## 5. Interpretation
CHRNB2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
