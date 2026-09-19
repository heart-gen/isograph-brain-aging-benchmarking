# SLC4A8 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.18 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.462  ·  missense o/e 0.62

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.18 | no | — | yes | — | risk allele G decreases SLC4A8 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Spinal_cord_cervical_c-1 | rs2241545 | G | -0.5919456481933594 | risk allele G decreases expression of SLC4A8 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS1, PABPC4, ENOX1, RBM41, A1CF, AKAP1, G3BP1, HNRNPA3, FXR1
- **All switched-motif RBPs (recurrence across regions):** AKAP1(4), AGO2(4), CELF4(4), CNOT4(4), CELF6(4), CELF5(4), DHX58(4), CPEB2(4), HNRNPA1L2(4), G3BP2(4), ESRP2(4), FXR1(4), ENOX1(4), EIF4B(4), HNRNPCL1(4), RBM24(4), RBM3(4), RBM41(4), RBM42(4), RBMS3(4)

## 5. Interpretation
SLC4A8 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
