# RUNDC3B — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.10 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.664  ·  missense o/e 0.94

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.10 | no | — | yes | — | risk allele C increases RUNDC3B expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs13233308 | C | 0.1919214129447937 | risk allele C increases expression of RUNDC3B |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, SNRPB2, RBMS3, PPRC1, HNRNPA3, CELF4, CELF5, CELF6, SNRNP70, YTHDC1, EIF4B, A1CF, CPEB2
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), ADAR(4), AGO1(4), AGO2(4), CELF4(4), CPEB2(4), CELF6(4), CELF5(4), CSTF2(4), DAZAP1(4), DDX19B(4), CPEB4(4), DDX58(4), EIF4A3(4), DHX9(4), DHX58(4), ZFP36L2(4), QKI(4), RALY(4), EIF4B(4)

## 5. Interpretation
RUNDC3B colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
