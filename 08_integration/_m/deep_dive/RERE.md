# RERE — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cortex)
- **Constraint:** LOEUF 0.352  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.02 | no | — | no | — | risk allele C increases RERE expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs3765971 | C | 0.3142753839492798 | risk allele C increases expression of RERE |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP2, ACO1, SAMD4A, RBM14, RBM42, A1CF, IGF2BP1, RBM6, G3BP1, RBMS3, PABPN1, RBM41, SRSF11, SNRPB2, RBM25
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), ADAR(3), AGO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CMTR1(3), CNOT4(3), CSTF2(3), CPEB2(3), DDX19B(3), DDX58(3), EIF4B(3), DHX58(3), DHX9(3), EIF4A3(3), ERI1(3), ENOX1(3)

## 5. Interpretation
RERE colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
