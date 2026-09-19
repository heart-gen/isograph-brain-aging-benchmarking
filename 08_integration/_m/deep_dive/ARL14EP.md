# ARL14EP — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Substantia_nigra)
- **Constraint:** LOEUF 0.823  ·  missense o/e 0.79

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Substantia_nigra | 0.03 | no | — | yes | — | risk allele A increases ARL14EP expression (gene-level; no intron) |
| scz | eQTL | Brain_Substantia_nigra | 0.03 | no | — | yes | — | risk allele A increases ARL14EP expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Substantia_nigra | rs1765142 | A | 0.28805291652679443 | risk allele A increases expression of ARL14EP |
| scz | Brain_Substantia_nigra | rs1765142 | A | 0.28805291652679443 | risk allele A increases expression of ARL14EP |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, RBM6, SNRPB2, RBMS3, G3BP2, ENOX1, IGF2BP1, CNOT4, HNRNPA3, CELF4, MATR3, PABPC3, RBM46
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), ADAR(6), AGO1(6), AGO2(6), CELF1(6), CELF2(6), CELF4(6), CMTR1(6), CNOT4(6), CPEB2(6), CPEB4(6), CSTF2(6), DDX19B(6), DDX58(6), DHX58(6), DHX9(6), EIF4A3(6), EIF4B(6), ENOX1(6)

## 5. Interpretation
ARL14EP colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
