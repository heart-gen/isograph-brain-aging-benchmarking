# MYO19 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.230  ·  missense o/e 1.01

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellum | 0.04 | no | — | no | — | risk allele A decreases MYO19 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellum | rs2285640 | A | -0.17561210691928864 | risk allele A decreases expression of MYO19 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM6, G3BP2, ENOX1, ZC3H10, IGF2BP1, HNRNPA3, CELF4, CELF5, CELF6, MATR3, SNRNP70, PABPC3, PABPC5, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), AGO1(4), CSTF2(4), AGO2(4), CELF4(4), CELF5(4), CPEB2(4), DHX58(4), DDX58(4), DDX19B(4), ZC3H10(4), ZRANB2(4), YBX2(4), EIF4A3(4), EIF4B(4), ENOX1(4), ESRP1(4), ESRP2(4), HNRNPAB(4), FXR2(4)

## 5. Interpretation
MYO19 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
