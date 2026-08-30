# WBP2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.769  ·  missense o/e 0.93

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | yes | — | risk allele T decreases WBP2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs74410877 | T | -0.33154040575027466 | risk allele T decreases expression of WBP2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** AGO1(5), AGO2(5), CELF5(5), CELF4(5), CSTF2(5), CPEB4(5), U2AF2(5), ZFP36L2(5), DDX19B(5), DHX58(5), EIF4A3(5), ELAVL3(5), ESRP1(5), HNRNPAB(5), HNRNPA3(5), HNRNPCL1(5), KHDRBS3(5), LIN28A(5), HNRNPM(5), HNRNPU(5)

## 5. Interpretation
WBP2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
