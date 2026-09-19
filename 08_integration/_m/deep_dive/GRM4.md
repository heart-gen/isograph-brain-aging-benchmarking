# GRM4 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.456  ·  missense o/e 0.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T increases GRM4 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs4713736 | T | 0.2643463909626007 | risk allele T increases expression of GRM4 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM6, SNRPB2, RBMS3, IGF2BP1, CNOT4, HNRNPA3, FXR1, MATR3, SNRNP70, PABPC3, RBM46, PABPC5, RBMS1, A1CF
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), AGO1(3), AGO2(3), AKAP1(3), CNOT4(3), CPEB1(3), CPEB4(3), DDX19B(3), DHX58(3), IGHMBP2(3), ELAVL3(3), ESRP2(3), FXR1(3), G3BP1(3), HNRNPA1L2(3), HNRNPDL(3), HNRNPM(3), IGF2BP2(3), HNRNPU(3), PABPC1(3)

## 5. Interpretation
GRM4 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
