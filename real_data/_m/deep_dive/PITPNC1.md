# PITPNC1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.418  ·  missense o/e 0.65

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Caudate_basal_ganglia | 0.02 | no | — | yes | — | risk allele C increases PITPNC1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs12103469 | C | 0.20539891719818115 | risk allele C increases expression of PITPNC1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, SART3, PABPC4, PUM2, RBMS3, EIF4A3, PABPC3, IFIH1, SRSF11, HNRNPCL1, CPEB2, DHX58, RALY, PHAX
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO1(4), AGO2(4), AKAP1(4), ANKHD1(4), CELF4(4), CELF5(4), CELF6(4), CMTR1(4), CPEB2(4), CSTF2(4), DDX58(4), DHX58(4), DHX9(4), EIF4A3(4), EIF4B(4), FXR2(4), G3BP1(4)

## 5. Interpretation
PITPNC1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
