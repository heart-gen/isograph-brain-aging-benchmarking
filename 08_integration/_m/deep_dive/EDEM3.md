# EDEM3 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.94 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.671  ·  missense o/e 0.87

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.94 | no | — | yes | — | risk allele G increases EDEM3 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs78444298 | G | 0.5194917321205139 | risk allele G increases expression of EDEM3 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, G3BP1, RBM6, SNRPB2, RBMS3, PPRC1, ENOX1, CNOT4, ZC3H10, IGF2BP1, HNRNPA3, CELF4, CELF5, CELF6, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** ESRP2(4), HNRNPA3(4), FXR2(4), RBM6(4), RBM8A(4), PUM2(4), A1CF(3), CELF4(3), ADAR(3), ENOX1(3), CELF5(3), HNRNPAB(3), ESRP1(3), DHX9(3), G3BP1(3), IGF2BP2(3), CELF6(3), DDX58(3), DHX58(3), YTHDC1(3)

## 5. Interpretation
EDEM3 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
