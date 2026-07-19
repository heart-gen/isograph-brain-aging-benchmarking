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
- **Switched *and* module-enriched (q<0.05) RBPs:** PPIE, KHDRBS1
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), ADAR(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CNOT4(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), DHX58(3), EIF4A3(3), EIF4B(3), ELAVL1(3), ELAVL3(3), ELAVL4(3), ENOX1(3)

## 5. Interpretation
DPYSL5 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
