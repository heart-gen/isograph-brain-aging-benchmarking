# WBP2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.769  ·  missense o/e 0.93

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | no | — | risk allele T decreases WBP2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs74410877 | T | -0.33154040575027466 | risk allele T decreases expression of WBP2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** AGO1(1), CELF4(1), CELF5(1), CPEB4(1), CSTF2(1), DDX19B(1), DHX58(1), EIF4A3(1), ELAVL3(1), ESRP1(1), HNRNPA3(1), HNRNPAB(1), HNRNPCL1(1), HNRNPLL(1), HNRNPM(1), HNRNPU(1), IGF2BP2(1), KHDRBS2(1), KHDRBS3(1), LIN28A(1)

## 5. Interpretation
WBP2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
