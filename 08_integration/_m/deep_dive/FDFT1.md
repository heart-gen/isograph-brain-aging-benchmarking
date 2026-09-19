# FDFT1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.07 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.701  ·  missense o/e 1.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Frontal_Cortex_BA9 | 0.07 | no | — | yes | — | risk allele G increases FDFT1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Frontal_Cortex_BA9 | rs2645429 | G | 0.2720873951911926 | risk allele G increases expression of FDFT1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, SNRPB2, RBMS3, CNOT4, CELF4, CELF5, PABPC3, RBM46, PABPC5, YTHDC1, EIF4B, A1CF, CPEB2, ZFP36L2, ZNF638
- **All switched-motif RBPs (recurrence across regions):** A1CF(7), AGO2(7), CELF4(7), CELF5(7), CSTF2(7), HNRNPLL(7), GRSF1(7), EIF4B(7), RBM28(7), SART3(7), SNRPB2(7), ZNF638(7), TRA2A(7), YTHDC1(7), SUPV3L1(7), SRSF10(7), RBFOX1(7), RBFOX2(7), RALY(7), PABPC5(7)

## 5. Interpretation
FDFT1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
