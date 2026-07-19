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
- **Switched *and* module-enriched (q<0.05) RBPs:** SUPV3L1, ZCRB1, RBM41, CNOT4
- **All switched-motif RBPs (recurrence across regions):** SUPV3L1(3), ZC3H10(3), RBM25(3), CELF6(2), ADAR(2), CELF4(2), CELF5(2), DHX58(2), DDX58(2), CNOT4(2), DHX9(2), IGF2BP1(2), IGF2BP2(2), G3BP1(2), ESRP1(2), SAMD4A(2), PABPN1(2), RBM28(2), RBM41(2), SNRPB2(2)

## 5. Interpretation
EDEM3 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
