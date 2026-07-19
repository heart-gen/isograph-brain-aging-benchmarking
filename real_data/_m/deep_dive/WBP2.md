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
- **All switched-motif RBPs (recurrence across regions):** AGO1(3), CPEB4(3), AGO2(3), DDX19B(3), CSTF2(3), ELAVL3(3), ESRP1(3), DHX58(3), EIF4A3(3), HNRNPLL(3), HNRNPCL1(3), HNRNPAB(3), HNRNPA3(3), SRSF10(3), SYNCRIP(3), YTHDC1(3), U2AF2(3), IGF2BP2(3), IGF2BP1(3), IGHMBP2(3)

## 5. Interpretation
WBP2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
