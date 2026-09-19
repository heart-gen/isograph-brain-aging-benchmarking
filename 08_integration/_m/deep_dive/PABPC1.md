# PABPC1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.09 (Brain_Hypothalamus)
- **Constraint:** LOEUF 0.119  ·  missense o/e 0.35

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Hypothalamus | 0.09 | no | — | no | — | risk allele C increases PABPC1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Hypothalamus | rs1693551 | C | 0.19282865524291992 | risk allele C increases expression of PABPC1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SAMD4A, RBM42, ZFP36L2, A1CF, IGF2BP1, G3BP1, PABPN1, SRSF11, SNRPB2, RBM25, PABPC4, DHX9, CSTF2, YTHDC1, ADAR
- **All switched-motif RBPs (recurrence across regions):** A1CF(7), AGO2(7), CSTF2(7), CPEB2(7), DDX58(7), DDX19B(7), G3BP1(7), EIF4A3(7), RBM25(7), KHDRBS2(7), KHDRBS3(7), HNRNPM(7), HNRNPLL(7), HNRNPCL1(7), SUPV3L1(7), TARDBP(7), U2AF2(7), YTHDC1(7), ZFP36L2(7), ZCRB1(7)

## 5. Interpretation
PABPC1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
