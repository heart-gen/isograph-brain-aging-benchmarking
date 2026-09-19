# SF3B1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cortex)
- **Constraint:** LOEUF 0.100  ·  missense o/e 0.28

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.02 | no | — | yes | — | risk allele C increases SF3B1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs12621129 | C | 0.2032228708267212 | risk allele C increases expression of SF3B1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM6, SNRPB2, PPRC1, ENOX1, ZC3H10, CNOT4, HNRNPA3, CELF4, CELF5, CELF6, MATR3, SNRNP70, PABPC3, PABPC5
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO1(4), AGO2(4), CELF4(4), CELF5(4), CELF6(4), CNOT4(4), CPEB2(4), DDX58(4), DDX19B(4), EIF4A3(4), DHX9(4), ZFP36L2(4), SNRPB2(4), ESRP1(4), ESRP2(4), G3BP1(4), FXR2(4)

## 5. Interpretation
SF3B1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
