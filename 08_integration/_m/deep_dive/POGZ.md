# POGZ — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 0.173  ·  missense o/e 0.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Putamen_basal_ganglia | 0.03 | no | — | yes | — | risk allele C decreases POGZ expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Putamen_basal_ganglia | rs6684085 | C | -0.17127901315689087 | risk allele C decreases expression of POGZ |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SNRPB2, ZC3H10, RBMS3, MSI1, RALY, RBM28, AKAP1, HNRNPCL1, CPEB2, SNRNP70, PABPC4, ZFP36L2, SAMD4A, PABPN1, SUPV3L1
- **All switched-motif RBPs (recurrence across regions):** AKAP1(6), CPEB2(6), MSI1(6), HNRNPA1L2(6), IGHMBP2(6), SNRPB2(6), SUPV3L1(6), RBMS3(6), RBM46(6), PUM2(6), PABPN1(6), PABPC4(6), SART3(6), ZFP36L2(2), HNRNPCL1(1), IGF2BP2(1), EIF4A3(1), DDX19B(1), CPEB1(1), HNRNPLL(1)

## 5. Interpretation
POGZ colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
