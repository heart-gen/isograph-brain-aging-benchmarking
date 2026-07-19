# CHRNB2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cortex)
- **Constraint:** LOEUF 1.116  ·  missense o/e 0.89

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.01 | no | — | no | — | risk allele C decreases CHRNB2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs11264227 | C | -0.12016735970973969 | risk allele C decreases expression of CHRNB2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PPIE, ELAVL4, KHDRBS1, HNRNPD, SYNCRIP, RNASEL, ELAVL3, OAS1, IGF2BP3, U2AF2, ELAVL2, PABPC1, KHDRBS3, RBMX, CELF2
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), ADAR(4), AGO1(4), AGO2(4), CELF2(4), CELF6(4), CPEB1(4), CPEB4(4), CSTF2(4), DAZAP1(4), DDX19B(4), DDX58(4), DHX58(4), DHX9(4), EIF4A3(4), ELAVL1(4), ELAVL2(4), ELAVL3(4), ELAVL4(4), ESRP1(4)

## 5. Interpretation
CHRNB2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
