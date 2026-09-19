# RCBTB1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.09 (Brain_Hypothalamus)
- **Constraint:** LOEUF 0.751  ·  missense o/e 0.82

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Hypothalamus | 0.09 | no | — | yes | — | risk allele G decreases RCBTB1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Hypothalamus | rs1925743 | G | -0.23535220324993134 | risk allele G decreases expression of RCBTB1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ACO1(2), ADAR(2), AGO1(2), CELF4(2), CELF5(2), CELF6(2), CMTR1(2), CPEB1(2), CPEB2(2), CPEB4(2), DAZAP1(2), DDX19B(2), DDX58(2), DHX58(2), DHX9(2), EIF4A3(2), ELAVL3(2), ENOX1(2), ESRP2(2), F2(2)

## 5. Interpretation
RCBTB1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
