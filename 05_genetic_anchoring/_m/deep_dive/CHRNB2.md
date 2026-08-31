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
- **Switched *and* module-enriched (q<0.05) RBPs:** ELAVL1, ELAVL4, KHDRBS1, RNASEL, PPIE, ZFP36, OAS1, HNRNPC, TIA1, IGF2BP3, HNRNPD, ELAVL2, PABPC1, RBMX, RC3H1
- **All switched-motif RBPs (recurrence across regions):** ADAR(5), AGO1(5), AGO2(5), CELF1(5), CELF2(5), CPEB1(5), CPEB4(5), CSTF2(5), DAZAP1(5), DDX19B(5), DDX58(5), DHX58(5), DHX9(5), EIF4A3(5), ELAVL1(5), ELAVL2(5), ELAVL3(5), ELAVL4(5), ENOX1(5), ESRP1(5)

## 5. Interpretation
CHRNB2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
