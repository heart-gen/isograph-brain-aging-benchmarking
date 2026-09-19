# MYO15A — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Amygdala)
- **Constraint:** LOEUF 0.889  ·  missense o/e 0.99

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Amygdala | 0.02 | no | — | no | — | risk allele T increases MYO15A expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Amygdala | rs854775 | T | 0.4060905873775482 | risk allele T increases expression of MYO15A |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBMS3, PABPC4, MATR3, CELF5, G3BP1, FXR1, AKAP1, YTHDC1, PTBP2
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), AGO2(2), AKAP1(2), CELF4(2), CELF5(2), CELF6(2), CPEB1(2), CPEB4(2), CSTF2(2), DHX58(2), EIF4B(2), ELAVL3(2), ENOX1(2), ESRP1(2), FXR1(2), G3BP1(2), HNRNPA1L2(2), HNRNPAB(2), HNRNPU(2)

## 5. Interpretation
MYO15A colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
