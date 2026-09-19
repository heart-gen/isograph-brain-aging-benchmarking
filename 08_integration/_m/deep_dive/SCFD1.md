# SCFD1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.697  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele G increases SCFD1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs229150 | G | 0.35154905915260315 | risk allele G increases expression of SCFD1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ACO1, SAMD4A, IGF2BP1, PABPN1, RBM25, CELF4, DHX9, CSTF2, YTHDC1, CELF5, SNRPB2, SART3, ADAR, IGHMBP2, MATR3
- **All switched-motif RBPs (recurrence across regions):** RBFOX1(5), LIN28A(5), A1CF(4), ACO1(4), AKAP1(4), CELF4(4), ADAR(4), AGO2(4), DHX9(4), PUM2(4), RBM25(4), NELFE(4), GRSF1(4), KHDRBS2(4), CELF5(4), CSTF2(4), SNRPA(4), TARDBP(4), ZFP36L2(4), ZRANB2(4)

## 5. Interpretation
SCFD1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
