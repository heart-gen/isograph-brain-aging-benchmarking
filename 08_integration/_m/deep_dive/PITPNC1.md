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
- **Switched *and* module-enriched (q<0.05) RBPs:** PUM1, SAMD4A, IFIH1, KHDRBS2, EIF4A3, PABPC4, G3BP1, RALY, RBMS3, HNRNPU, SF1, DHX58, SRSF11, SYNCRIP
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), ADAR(5), AGO1(5), AGO2(5), AKAP1(5), ANKHD1(5), CELF4(5), CELF5(5), CMTR1(5), CNOT4(5), CPEB2(5), CSTF2(5), DAZAP1(5), DDX19B(5), DDX58(5), DHX58(5), DHX9(5), EIF4A3(5), EIF4B(5)

## 5. Interpretation
PITPNC1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
