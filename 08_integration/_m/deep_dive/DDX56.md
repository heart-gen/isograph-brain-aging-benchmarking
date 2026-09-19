# DDX56 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.886  ·  missense o/e 1.00

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele C decreases DDX56 expression (gene-level; no intron) |
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele C decreases DDX56 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs217364 | C | -0.16198846697807312 | risk allele C decreases expression of DDX56 |
| scz | Brain_Cerebellar_Hemisphere | rs217364 | C | -0.16198846697807312 | risk allele C decreases expression of DDX56 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SRSF11, RBM25, CSTF2, SRSF10, PABPC1, SART3, IGHMBP2, CPEB4, MATR3, SNRNP70, HNRNPA3, ZNF638, ZRANB2, EIF4B, GRSF1
- **All switched-motif RBPs (recurrence across regions):** AGO1(5), AGO2(5), AKAP1(5), CPEB2(5), EIF4A3(5), DHX58(5), DDX19B(5), ELAVL3(5), HNRNPD(5), HNRNPU(5), HNRNPAB(5), HNRNPA1L2(5), HNRNPA3(5), HNRNPA0(5), ESRP1(5), SART3(5), SF1(5), TARDBP(5), ZRANB2(5), ZC3H10(5)

## 5. Interpretation
DDX56 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
