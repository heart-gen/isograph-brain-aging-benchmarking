# FBXL15 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.048  ·  missense o/e 1.23

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | no | — | risk allele C increases FBXL15 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs117842157 | C | 0.745553195476532 | risk allele C increases expression of FBXL15 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPC, IGF2BP3, KHDRBS1, ELAVL4, DAZAP1, HNRNPU, HNRNPK, TIAL1, SYNCRIP, PABPC1, ELAVL1, ZFP36, NUDT21, QKI
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), AGO1(3), AGO2(3), CELF2(3), CELF4(3), CNOT4(3), CPEB1(3), CPEB2(3), CPEB4(3), DAZAP1(3), DDX19B(3), DHX58(3), EIF4A3(3), EIF4B(3), ELAVL1(3), ELAVL3(3), ELAVL4(3), FXR2(3), GRSF1(3), HNRNPA3(3)

## 5. Interpretation
FBXL15 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
