# CR1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.17 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.774  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.17 | no | — | yes | — | risk allele A increases CR1 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** U2AF2, SYNCRIP, HNRNPU, QKI, KHDRBS3, RBM4, DDX19B, HNRNPDL
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), AGO2(3), CELF4(3), CELF5(3), CPEB2(3), CELF6(3), CSTF2(3), DDX19B(3), FXR2(3), DHX58(3), EIF4A3(3), EIF4B(3), GRSF1(3), G3BP1(3), HNRNPCL1(3), HNRNPA3(3), ZFP36L2(3), U2AF2(3), KHDRBS2(3), IGHMBP2(3)

## 5. Interpretation
CR1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
