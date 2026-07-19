# STK39 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.04 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.381  ·  missense o/e 0.72

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.04 | no | — | yes | — | risk allele T decreases STK39 expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, HNRNPCL1, PABPC5, A1CF, MATR3, RBMS3, ZFP36L2, G3BP2, PABPC4, RBMS1, HNRNPA1L2, RBM6, ESRP1, AKAP1, ADAR
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO2(4), AKAP1(4), CPEB1(4), DDX58(4), DHX58(4), DHX9(4), EIF4A3(4), ESRP1(4), G3BP1(4), G3BP2(4), HNRNPA1L2(4), HNRNPA3(4), HNRNPAB(4), HNRNPCL1(4), HNRNPM(4), IFIH1(4), MATR3(4)

## 5. Interpretation
STK39 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
