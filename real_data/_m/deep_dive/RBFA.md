# RBFA — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.977  ·  missense o/e 1.03

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele C increases usage of junction chr18:80042219-80044212(+) (ENST00000306735.10,E |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr18:80048193-80069798(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele C decreases usage of junction chr18:80067308-80069798(+) (unmapped transcript; |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** IGF2BP1
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ACO1(1), AGO1(1), AKAP1(1), ANKHD1(1), CELF4(1), CELF5(1), CELF6(1), CMTR1(1), CPEB1(1), CPEB2(1), CPEB4(1), DDX19B(1), DHX58(1), EIF4A3(1), EIF4B(1), ELAVL3(1), ENOX1(1), ESRP1(1), ESRP2(1)

## 5. Interpretation
RBFA has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
