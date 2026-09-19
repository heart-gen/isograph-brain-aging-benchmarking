# TYK2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** LBD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.622  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.03 | no | — | no | — | risk allele A decreases TYK2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| lbd | Brain_Spinal_cord_cervical_c-1 | rs8101195 | A | -0.19877827167510986 | risk allele A decreases expression of TYK2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBMS3, ENOX1, CNOT4, IGF2BP1, HNRNPA3, CELF4, CELF5, CELF6, MATR3, PABPC5, YTHDC1, A1CF, CPEB2, ZFP36L2
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AKAP1(3), CELF6(3), CNOT4(3), CPEB4(3), CPEB1(3), CPEB2(3), HNRNPA0(3), FXR2(3), ESRP2(3), EIF4A3(3), DDX19B(3), KHDRBS1(3), KHDRBS2(3), IGF2BP1(3), HNRNPAB(3), HNRNPA3(3), GRSF1(3), SYNCRIP(3)

## 5. Interpretation
TYK2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
