# ARHGAP44 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Amygdala)
- **Constraint:** LOEUF 0.587  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Amygdala | 0.02 | no | — | yes | — | risk allele A decreases ARHGAP44 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Amygdala | rs12936683 | A | -0.2662936747074127 | risk allele A decreases expression of ARHGAP44 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SART3, RBM8A, IGF2BP2, PPRC1, SUPV3L1, SNRPB2, HNRNPAB, AGO1, ZFP36L2, AKAP1, AGO2, PABPC4, A1CF, EIF4A3, ENOX1
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO2(3), AKAP1(3), EIF4B(3), CELF6(3), CNOT4(3), CPEB2(3), CPEB4(3), CSTF2(3), DHX58(3), EIF4A3(3), G3BP1(3), ENOX1(3), ESRP1(3), ESRP2(3), HNRNPA0(3), G3BP2(3), HNRNPA3(3), HNRNPAB(3)

## 5. Interpretation
ARHGAP44 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
