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
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM6, ZC3H10, RBM14, FXR2, RBMS3, G3BP1, SNRPB2, A1CF, HNRNPA3, IGF2BP1, RBM28, ENOX1, ZCRB1, PABPN1, RBM8A
- **All switched-motif RBPs (recurrence across regions):** FXR2(6), HNRNPA3(6), RBM8A(6), A1CF(3), ADAR(3), CELF4(3), DDX58(3), DHX58(3), CELF6(3), CELF5(3), ENOX1(3), DHX9(3), ESRP2(3), ESRP1(3), G3BP1(3), HNRNPAB(3), IGF2BP1(3), CNOT4(3), IGF2BP2(3), NELFE(3)

## 5. Interpretation
EDEM3 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
