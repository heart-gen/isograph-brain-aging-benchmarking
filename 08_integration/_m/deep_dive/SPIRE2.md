# SPIRE2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.346  ·  missense o/e 1.28

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele A decreases SPIRE2 expression (gene-level; no intron) |
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele A decreases SPIRE2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs11076631 | A | -0.3170987665653229 | risk allele A decreases expression of SPIRE2 |
| scz | Brain_Cerebellar_Hemisphere | rs11076631 | A | -0.3170987665653229 | risk allele A decreases expression of SPIRE2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** IGF2BP3, HNRNPD, ELAVL3, KHDRBS1, AGO1, RBMY1A1, NOVA2, PABPC1, PUM1, YTHDC1, HNRNPDL
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), AGO1(2), AGO2(2), AKAP1(2), CELF1(2), CELF2(2), CELF4(2), CELF5(2), CELF6(2), CPEB1(2), EIF4A3(2), CPEB2(2), CSTF2(2), DDX19B(2), EIF4B(2), DHX58(2), ELAVL4(2), ELAVL3(2), HNRNPDL(2)

## 5. Interpretation
SPIRE2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
