# ANKRD54 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.885  ·  missense o/e 0.75

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.02 | no | — | yes | — | risk allele A decreases ANKRD54 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Spinal_cord_cervical_c-1 | rs74622487 | A | -0.621799647808075 | risk allele A decreases expression of ANKRD54 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, HNRNPA0, HNRNPC, ELAVL3, KHDRBS1, HNRNPD, ELAVL4, IGF2BP3, CPEB1, RNASEL, SRSF10, PABPC1, CPEB4, ZFP36, TIAL1
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), CELF1(3), ELAVL3(3), CELF2(3), CELF6(3), CPEB1(3), CPEB4(3), CSTF2(3), EIF4A3(3), DDX19B(3), ERI1(3), ELAVL4(3), EIF4B(3), ELAVL1(3), ESRP2(3), HNRNPA0(3), GRSF1(3), FXR1(3), TRA2A(3)

## 5. Interpretation
ANKRD54 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
