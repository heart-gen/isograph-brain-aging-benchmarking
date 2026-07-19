# FDFT1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.07 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.701  ·  missense o/e 1.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Frontal_Cortex_BA9 | 0.07 | no | — | yes | — | risk allele G increases FDFT1 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB1
- **All switched-motif RBPs (recurrence across regions):** AKAP1(3), CELF5(3), PABPC5(3), HNRNPA3(3), SRSF10(3), SNRPB2(3), RBM28(3), RBFOX2(3), YTHDC1(3), SUPV3L1(3), ZNF638(3), DHX58(1), NUDT21(1), NELFE(1), CPEB1(1), RALY(1)

## 5. Interpretation
FDFT1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
