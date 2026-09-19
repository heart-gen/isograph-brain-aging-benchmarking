# CNNM2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.287  ·  missense o/e 0.59

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | no | — | risk allele T decreases CNNM2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs11191580 | T | -0.4082171618938446 | risk allele T decreases expression of CNNM2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PUM1, SF1, SRSF6, TARDBP, FMR1, NONO, CPEB4, EIF4A3, DDX19B
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), AGO2(2), AKAP1(2), ANKHD1(2), CELF4(2), CELF5(2), CELF6(2), CNOT4(2), CPEB1(2), CPEB2(2), CPEB4(2), CSTF2(2), DDX19B(2), DHX58(2), EIF4A3(2), EIF4B(2), ENOX1(2), ESRP1(2), ESRP2(2)

## 5. Interpretation
CNNM2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
