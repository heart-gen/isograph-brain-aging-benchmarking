# PKD1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.360  ·  missense o/e 1.22

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele G decreases usage of junction chr16:2100566-2102102(-) (unmapped transcript; n |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele G decreases usage of junction chr16:2105662-2105865(-) (ENST00000564865.5; not |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele G increases usage of junction chr16:2105935-2106091(-) (unmapped transcript; n |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele G decreases usage of junction chr16:2106024-2106091(-) (ENST00000262304.9,ENST |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele G decreases usage of junction chr16:2112399-2112788(-) (unmapped transcript; n |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM14, RBM41, PPRC1, ZC3H10, ENOX1, RBMS3, CELF6, RALY, CNOT4, ZCRB1, MATR3, PABPC3, A1CF, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), AGO1(6), AGO2(6), AKAP1(6), CELF4(6), CELF5(6), CELF6(6), CNOT4(6), CPEB2(6), DAZAP1(6), CPEB4(6), DHX58(6), DDX19B(6), ELAVL3(6), ELAVL4(6), EIF4A3(6), ELAVL2(6), ENOX1(6), ESRP1(6)

## 5. Interpretation
PKD1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
