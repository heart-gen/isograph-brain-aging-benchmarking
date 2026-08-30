# ADAM10 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.05 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.276  ·  missense o/e 0.65

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.05 | no | — | yes | — | risk allele T increases ADAM10 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Spinal_cord_cervical_c-1 | rs1427281 | T | 0.3229811489582062 | risk allele T increases expression of ADAM10 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM41, SNRPB2, PPRC1, RBM6, RBMS3, CELF6, MSI1, RALY, CNOT4, ZCRB1, RBMS1, RBM28, A1CF, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO2(4), AKAP1(4), CELF4(4), CELF5(4), G3BP2(4), CELF6(4), CNOT4(4), CPEB2(4), CSTF2(4), DDX19B(4), DHX9(4), EIF4A3(4), ESRP2(4), ESRP1(4), FXR1(4), HNRNPU(4), HNRNPLL(4)

## 5. Interpretation
ADAM10 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
