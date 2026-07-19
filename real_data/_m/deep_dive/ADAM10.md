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
- **Switched *and* module-enriched (q<0.05) RBPs:** ZCRB1, RALY, HNRNPCL1, CPEB2, PABPC5, A1CF, MATR3, YBX2, CELF6, RBMS3, ZFP36L2, G3BP2, PABPC4, IGF2BP1, FXR1
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), ADAR(3), CELF5(3), CELF4(3), CELF6(3), CPEB2(3), DDX19B(3), CSTF2(3), EIF4A3(3), EIF4B(3), DHX58(3), DHX9(3), ESRP1(3), G3BP2(3), FXR1(3), ESRP2(3), IGF2BP2(3), KHDRBS2(3), KHDRBS3(3)

## 5. Interpretation
ADAM10 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
