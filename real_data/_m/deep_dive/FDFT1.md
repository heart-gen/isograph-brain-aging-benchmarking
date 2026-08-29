# FDFT1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.07 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.701  ·  missense o/e 1.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Frontal_Cortex_BA9 | 0.07 | no | — | yes | — | risk allele G increases FDFT1 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB1, PUM1, NONO
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), AGO2(3), CELF4(3), CELF5(3), GRSF1(3), CSTF2(3), HNRNPLL(3), RBM28(3), SNRPA(3), RBMS3(3), SRSF10(3), YTHDC1(3), ZNF638(3), SUPV3L1(3), TRA2A(3), SNRPB2(3), RBFOX1(3), RBFOX2(3), PABPC1(3), PABPC5(3)

## 5. Interpretation
FDFT1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
