# RBM26 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.382  ·  missense o/e 0.79

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | — | no | — | risk allele A increases RBM26 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs9318627 | A | 0.2170882523059845 | risk allele A increases expression of RBM26 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, RBM6, RBMS3, PPRC1, IGF2BP1, HNRNPA3, CELF4, CELF5, FXR1, CELF6, PABPC3, RBM46, PABPC5
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), CELF5(4), CELF4(4), CELF6(4), DHX58(4), CSTF2(4), CPEB2(4), SART3(4), YTHDC1(4), SUPV3L1(4), EIF4B(4), ESRP1(4), FXR2(4), FXR1(4), ESRP2(4), HNRNPAB(4), G3BP1(4), PABPC4(4), PABPC3(4), KHDRBS2(4)

## 5. Interpretation
RBM26 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
