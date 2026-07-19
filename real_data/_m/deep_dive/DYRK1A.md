# DYRK1A — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.172  ·  missense o/e 0.64

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele T increases DYRK1A expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZCRB1, RBM41, RALY, HNRNPCL1, CPEB2, PABPC5, A1CF, RBM42, MATR3, PABPC3, CELF6, RBMS3, ZFP36L2, G3BP2, PABPC4
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CPEB2(3), DDX19B(3), EIF4A3(3), EIF4B(3), ENOX1(3), ESRP1(3), ESRP2(3), FXR2(3), G3BP2(3), HNRNPA0(3), HNRNPA1L2(3), HNRNPCL1(3), KHDRBS2(3), KHDRBS3(3), MATR3(3)

## 5. Interpretation
DYRK1A colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
