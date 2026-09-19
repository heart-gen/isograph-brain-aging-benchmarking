# PCCB — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cortex)
- **Constraint:** LOEUF 0.843  ·  missense o/e 1.11

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.01 | no | — | yes | — | risk allele C decreases PCCB expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs66691851 | C | -0.4196547865867615 | risk allele C decreases expression of PCCB |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM6, SNRPB2, PPRC1, G3BP2, IGF2BP1, HNRNPA3, CELF6, MATR3, PABPC3, RBM46, PABPC5, YTHDC1, A1CF, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), AGO2(3), CELF6(3), CPEB2(3), ELAVL3(3), FXR2(3), ESRP1(3), G3BP1(3), G3BP2(3), HNRNPA0(3), GRSF1(3), HNRNPU(3), HNRNPDL(3), HNRNPAB(3), HNRNPD(3), KHDRBS3(3), KHDRBS2(3), IGHMBP2(3), IGF2BP1(3), SART3(3)

## 5. Interpretation
PCCB colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
