# DYRK1A — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.172  ·  missense o/e 0.64

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele T increases DYRK1A expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM14, RBM41, SNRPB2, PPRC1, RBM6, ZC3H10, ENOX1, RBMS3, CELF6, RALY, RBMS1, ZCRB1, MATR3, PABPC3
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), AGO2(5), ANKHD1(5), CELF4(5), CELF5(5), CELF6(5), CPEB2(5), DDX19B(5), EIF4A3(5), EIF4B(5), ESRP1(5), ENOX1(5), ESRP2(5), FXR2(5), G3BP2(5), G3BP1(5), KHDRBS3(5), MATR3(5), GRSF1(5), HNRNPA0(5)

## 5. Interpretation
DYRK1A colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
