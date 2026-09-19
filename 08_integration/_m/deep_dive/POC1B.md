# POC1B — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Anterior_cingulate_cortex_BA24)
- **Constraint:** LOEUF n/a  ·  missense o/e n/a

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.03 | no | — | no | — | risk allele T increases POC1B expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs1492914 | T | 0.1853104829788208 | risk allele T increases expression of POC1B |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM42, PUM1, RBM6, G3BP1, SAMD4A, SNRPB2, CNOT4, RBM41, G3BP2, PPRC1, DDX58, TARDBP, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CNOT4(3), CPEB1(3), CPEB2(3), CSTF2(3), DDX19B(3), DDX58(3), DHX58(3), EIF4A3(3), ENOX1(3), ESRP1(3), ESRP2(3), FXR2(3), G3BP1(3)

## 5. Interpretation
POC1B colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
