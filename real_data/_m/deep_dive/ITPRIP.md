# ITPRIP — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.06 (Brain_Cortex)
- **Constraint:** LOEUF 1.140  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.06 | no | — | yes | — | risk allele G decreases ITPRIP expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs17883150 | G | -0.2046354115009308 | risk allele G decreases expression of ITPRIP |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CNOT4(3), CPEB1(3), CPEB2(3), CPEB4(3), DHX58(3), EIF4A3(3), ERI1(3), FXR2(3), G3BP2(3), HNRNPA0(3), HNRNPA1L2(3), HNRNPA3(3), HNRNPAB(3), HNRNPD(3), HNRNPLL(3)

## 5. Interpretation
ITPRIP colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
