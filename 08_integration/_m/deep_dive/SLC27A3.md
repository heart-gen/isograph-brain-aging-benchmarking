# SLC27A3 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.262  ·  missense o/e 1.10

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | — | yes | — | risk allele A decreases SLC27A3 expression (gene-level; no intron) |
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | — | no | — | risk allele A decreases SLC27A3 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs4540690 | A | -0.2447323203086853 | risk allele A decreases expression of SLC27A3 |
| scz | Brain_Cerebellar_Hemisphere | rs4540690 | A | -0.2447323203086853 | risk allele A decreases expression of SLC27A3 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM6, SNRPB2, RBMS3, PPRC1, G3BP2, ENOX1, ZC3H10, HNRNPA3, CELF6, MATR3, PABPC3, PABPC5, YTHDC1, RBMS1, A1CF
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), CELF6(6), CPEB1(6), CPEB2(6), CPEB4(6), CSTF2(6), DDX58(6), DHX58(6), ELAVL3(6), ENOX1(6), ESRP2(6), FMR1(6), G3BP2(6), HNRNPA0(6), HNRNPA3(6), HNRNPAB(6), HNRNPC(6), HNRNPD(6), HNRNPDL(6)

## 5. Interpretation
SLC27A3 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
