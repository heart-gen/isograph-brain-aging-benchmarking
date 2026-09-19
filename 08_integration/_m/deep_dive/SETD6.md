# SETD6 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cortex)
- **Constraint:** LOEUF 1.617  ·  missense o/e 1.22

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.01 | no | — | yes | — | risk allele T decreases SETD6 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs4784053 | T | -0.2903788685798645 | risk allele T decreases expression of SETD6 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM6, PPRC1, G3BP2, IGF2BP1, CNOT4, HNRNPA3, CELF4, CELF5, CELF6, RBM46, PABPC5, YTHDC1, EIF4B, A1CF
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), AGO1(5), AGO2(5), CELF4(5), CELF5(5), CELF6(5), CNOT4(5), CPEB2(5), CPEB4(5), ZCRB1(5), CSTF2(5), DDX19B(5), DHX58(5), EIF4B(5), EIF4A3(5), G3BP2(5), ESRP1(5), HNRNPCL1(5), HNRNPA3(5), HNRNPAB(5)

## 5. Interpretation
SETD6 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
