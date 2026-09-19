# OLA1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cortex)
- **Constraint:** LOEUF 0.434  ·  missense o/e 0.75

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.01 | no | — | yes | — | risk allele A decreases OLA1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs35084581 | A | -0.2617471516132355 | risk allele A decreases expression of OLA1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PUM1, SNRPB2, RBMS1, PPRC1, NUDT21, RBM5, ERI1, ZRANB2, PABPC4, SRSF6, G3BP2, HNRNPDL, TARDBP, AKAP1, NOVA2
- **All switched-motif RBPs (recurrence across regions):** ACO1(6), AGO1(6), AKAP1(6), AGO2(6), CELF6(6), CPEB2(6), DDX58(6), CPEB4(6), CSTF2(6), DDX19B(6), ESRP2(6), ESRP1(6), ERI1(6), EIF4A3(6), YTHDC1(6), ZFP36L2(6), SNRPA(6), SNRPB2(6), FXR1(6), FUS(6)

## 5. Interpretation
OLA1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
