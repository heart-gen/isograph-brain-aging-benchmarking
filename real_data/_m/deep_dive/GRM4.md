# GRM4 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.456  ·  missense o/e 0.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele T increases GRM4 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs4713736 | T | 0.2643463909626007 | risk allele T increases expression of GRM4 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** DHX58
- **All switched-motif RBPs (recurrence across regions):** AGO2(3), CNOT4(3), CPEB4(3), CPEB1(3), HNRNPA1L2(3), DHX58(3), SART3(3), SNRNP70(3), G3BP1(3), FXR2(3), PABPC4(3), SAMD4A(3), RBFOX2(3), PABPC5(3), SUPV3L1(3), AGO1(1), AKAP1(1), HNRNPA3(1), HNRNPA0(1), ESRP2(1)

## 5. Interpretation
GRM4 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
