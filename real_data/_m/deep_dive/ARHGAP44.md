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
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB4, RBM14, RBMS3, SUPV3L1, CELF4, RBM24, CSTF2
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), AGO2(4), AKAP1(4), CELF5(4), CELF4(4), CELF6(4), CPEB2(4), DHX58(4), CPEB4(4), CSTF2(4), DAZAP1(4), EIF4A3(4), DDX19B(4), ELAVL3(4), EIF4B(4), ZCRB1(4), RALY(4), ENOX1(4), ESRP1(4)

## 5. Interpretation
ARHGAP44 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
