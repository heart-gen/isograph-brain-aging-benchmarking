# SH3GL2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.05 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.937  ·  missense o/e 1.01

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Frontal_Cortex_BA9 | 0.05 | no | — | no | — | risk allele T increases SH3GL2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Frontal_Cortex_BA9 | rs10756905 | T | 0.1845596730709076 | risk allele T increases expression of SH3GL2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM6, ZC3H10, CELF4, CELF5, FXR1, PABPC3, RBMS1, ZFP36L2, IGF2BP2, PABPC4, RALY, NELFE, HNRNPM, HNRNPLL
- **All switched-motif RBPs (recurrence across regions):** ACO1(7), AGO1(7), CELF4(7), CELF5(7), CMTR1(7), CPEB1(7), CSTF2(7), DAZAP1(7), DDX19B(7), EIF4A3(7), ELAVL1(7), ELAVL2(7), ELAVL3(7), ELAVL4(7), FMR1(7), FXR1(7), GRSF1(7), HNRNPA0(7), HNRNPA1L2(7), HNRNPA2B1(7)

## 5. Interpretation
SH3GL2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
