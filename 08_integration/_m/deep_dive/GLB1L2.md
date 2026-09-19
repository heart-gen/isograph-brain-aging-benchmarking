# GLB1L2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.05 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.021  ·  missense o/e 0.91

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Putamen_basal_ganglia | 0.05 | no | — | yes | — | risk allele C decreases GLB1L2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Putamen_basal_ganglia | rs7112912 | C | -0.22864803671836853 | risk allele C decreases expression of GLB1L2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PABPN1, PUM2, CELF6, RBM25, SUPV3L1, IGF2BP2, KHDRBS1, PABPC1, DHX58, ACO1, PABPC3
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ACO1(1), ADAR(1), AGO1(1), AGO2(1), CELF6(1), CMTR1(1), CNOT4(1), CPEB1(1), CPEB2(1), CPEB4(1), CSTF2(1), DAZAP1(1), DDX19B(1), DHX58(1), DHX9(1), EIF4A3(1), ELAVL3(1), ESRP1(1), FXR1(1)

## 5. Interpretation
GLB1L2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
