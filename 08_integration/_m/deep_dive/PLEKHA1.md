# PLEKHA1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.916  ·  missense o/e 1.00

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | — | yes | — | risk allele T increases PLEKHA1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellar_Hemisphere | rs56029211 | T | 0.22710636258125305 | risk allele T increases expression of PLEKHA1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, RBM6, SNRPB2, RBMS3, PPRC1, IGF2BP1, ZC3H10, CNOT4, HNRNPA3, CELF4, CELF5, CELF6, MATR3
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), AGO2(5), AKAP1(5), CELF5(5), CNOT4(5), ESRP1(5), DHX58(5), PABPC4(5), MATR3(5), HNRNPM(5), IGF2BP2(5), IGF2BP1(5), IGHMBP2(5), HNRNPLL(5), FXR2(5), ESRP2(5), SNRPB2(5), SRSF4(5), SUPV3L1(5), YTHDC1(5)

## 5. Interpretation
PLEKHA1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
