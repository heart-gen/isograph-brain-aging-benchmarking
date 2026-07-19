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
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB2, PABPC3, RBMS3, PABPC4, PABPN1, HNRNPA1L2, RBM46, SART3, CELF4
- **All switched-motif RBPs (recurrence across regions):** CELF4(3), CPEB2(3), HNRNPA1L2(3), PABPC3(3), PABPC4(3), PABPN1(3), RBM46(3), RBMS3(3), SAMD4A(3), SART3(3), SNRPB2(3), SUPV3L1(3)

## 5. Interpretation
POGZ colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
