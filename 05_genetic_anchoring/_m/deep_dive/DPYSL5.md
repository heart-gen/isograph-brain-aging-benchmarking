# DPYSL5 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Substantia_nigra)
- **Constraint:** LOEUF 0.268  ·  missense o/e 0.62

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Substantia_nigra | 0.01 | no | — | yes | — | risk allele G increases DPYSL5 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Substantia_nigra | rs2082610 | G | 0.21863797307014465 | risk allele G increases expression of DPYSL5 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ELAVL4, ELAVL1
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), ADAR(4), AGO1(4), AGO2(4), AKAP1(4), CELF4(4), CELF5(4), CELF6(4), CNOT4(4), CPEB1(4), CPEB2(4), CPEB4(4), CSTF2(4), DAZAP1(4), DDX19B(4), DHX58(4), EIF4A3(4), EIF4B(4), ELAVL1(4), ELAVL3(4)

## 5. Interpretation
DPYSL5 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
