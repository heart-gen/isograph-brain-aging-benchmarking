# CR1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.17 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.774  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.17 | no | — | yes | — | risk allele A increases CR1 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** QKI, AGO2, PHAX
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), CELF4(3), AKAP1(3), ZNF638(3), SAMD4A(3), CELF5(3), CPEB2(3), EIF4A3(3), DHX58(3), EIF4B(3), G3BP1(3), HNRNPCL1(3), MATR3(3), HNRNPA3(3), SART3(3), QKI(3), RALY(3), PHAX(3), PABPC4(3), RBMS1(3)

## 5. Interpretation
CR1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
