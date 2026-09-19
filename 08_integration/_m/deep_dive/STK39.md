# STK39 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.04 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.381  ·  missense o/e 0.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.04 | no | — | no | — | risk allele T decreases STK39 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Nucleus_accumbens_basal_ganglia | rs1474055 | T | -0.21550074219703674 | risk allele T decreases expression of STK39 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ACO1, SAMD4A, RBM14, ZFP36L2, A1CF, RBM6, G3BP1, SRSF11, RBM25, PABPC4, DHX9, CSTF2, ADAR, SART3, IGHMBP2
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), ADAR(6), AGO2(6), AKAP1(6), CELF2(6), CPEB1(6), CSTF2(6), DDX19B(6), DDX58(6), DHX58(6), DHX9(6), EIF4A3(6), ESRP1(6), G3BP1(6), GRSF1(6), HNRNPA1L2(6), HNRNPA3(6), HNRNPAB(6), HNRNPCL1(6)

## 5. Interpretation
STK39 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
