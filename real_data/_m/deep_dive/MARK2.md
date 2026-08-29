# MARK2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Cortex)
- **Constraint:** LOEUF 0.201  ·  missense o/e 0.57

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.03 | no | — | yes | — | risk allele T increases MARK2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs7121067 | T | 0.1440276950597763 | risk allele T increases expression of MARK2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, SNRPB2, RBMS3, A1CF, AKAP1, IGF2BP1, ZC3H10, ESRP1, RBM42, G3BP2, CNOT4, ACO1, CPEB2, PABPN1, HNRNPLL
- **All switched-motif RBPs (recurrence across regions):** A1CF(7), ACO1(7), AGO2(7), AKAP1(7), CELF4(7), CELF5(7), CELF6(7), CNOT4(7), CPEB1(7), CPEB2(7), CSTF2(7), DDX19B(7), DHX58(7), EIF4A3(7), EIF4B(7), ENOX1(7), G3BP1(7), HNRNPCL1(7), HNRNPU(7), HNRNPM(7)

## 5. Interpretation
MARK2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
