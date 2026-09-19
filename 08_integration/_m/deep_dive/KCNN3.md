# KCNN3 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Anterior_cingulate_cortex_BA24)
- **Constraint:** LOEUF 0.379  ·  missense o/e 0.60

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.03 | no | — | no | — | risk allele A decreases KCNN3 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs7531728 | A | -0.17171010375022888 | risk allele A decreases expression of KCNN3 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPC, HNRNPA0, KHDRBS1, ELAVL4, HNRNPU, TIAL1, SYNCRIP, QKI, U2AF2, ELAVL3, PABPC5
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), DDX58(3), DHX58(3), EIF4A3(3), ELAVL3(3), ELAVL4(3), ENOX1(3), ESRP1(3)

## 5. Interpretation
KCNN3 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
