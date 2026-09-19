# GALNT6 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.18 (Brain_Cortex)
- **Constraint:** LOEUF 1.024  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Cortex | 0.18 | no | — | yes | — | risk allele G decreases GALNT6 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cortex | rs2241545 | G | -0.5057262182235718 | risk allele G decreases expression of GALNT6 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPC, HNRNPA0, KHDRBS1, DAZAP1, HNRNPU, HNRNPK, TIAL1, SYNCRIP, PABPC1, NUDT21, QKI, U2AF2, ELAVL3
- **All switched-motif RBPs (recurrence across regions):** ACO1(2), AGO1(2), AGO2(2), AKAP1(2), CELF6(2), CNOT4(2), CMTR1(2), CPEB1(2), CPEB2(2), CSTF2(2), CPEB4(2), DHX58(2), EIF4A3(2), DAZAP1(2), DDX19B(2), ELAVL3(2), ESRP1(2), ERI1(2), ENOX1(2), ZNF638(2)

## 5. Interpretation
GALNT6 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
