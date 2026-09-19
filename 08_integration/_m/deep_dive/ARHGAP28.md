# ARHGAP28 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.07 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.773  ·  missense o/e 0.91

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.07 | no | — | yes | — | risk allele T decreases ARHGAP28 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs77530930 | T | -0.48435094952583313 | risk allele T decreases expression of ARHGAP28 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ADAR(1), AGO2(1), AKAP1(1), CELF4(1), CELF5(1), CELF6(1), CNOT4(1), CPEB1(1), CPEB2(1), CPEB4(1), CSTF2(1), DDX19B(1), DDX58(1), DHX58(1), DHX9(1), EIF4A3(1), ELAVL3(1), ESRP1(1), ESRP2(1)

## 5. Interpretation
ARHGAP28 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
