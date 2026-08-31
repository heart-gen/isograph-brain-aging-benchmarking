# STK39 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.04 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.381  ·  missense o/e 0.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.04 | no | — | yes | — | risk allele T decreases STK39 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM14, RBM6, ZC3H10, RBMS1, A1CF, AKAP1, ESRP1, HNRNPA3, HNRNPCL1, ADAR, DHX9, PABPC4, PABPC5, ZFP36L2
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), ADAR(5), AGO2(5), AKAP1(5), CELF2(5), CPEB1(5), CSTF2(5), DDX19B(5), DDX58(5), DHX58(5), DHX9(5), EIF4A3(5), ESRP1(5), G3BP1(5), GRSF1(5), HNRNPA1L2(5), HNRNPA3(5), HNRNPAB(5), HNRNPCL1(5)

## 5. Interpretation
STK39 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
