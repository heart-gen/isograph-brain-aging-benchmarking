# SLC36A1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** LBD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.737  ·  missense o/e 0.82

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | sQTL | Brain_Cerebellum | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr5:151447813-151458788(+) (ENST00000243389.8,E |
| lbd | sQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele G decreases usage of junction chr5:151452216-151458788(+) (unmapped transcript |
| lbd | sQTL | Brain_Cerebellum | 0.01 | no | no | no | — | risk allele G increases usage of junction chr5:151452491-151458788(+) (ENST00000517945.5,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| lbd | Brain_Cerebellum | rs6875844 | G | -0.9302634000778198 | risk allele G decreases intron usage of chr5:151447813-151458788(+) |
| lbd | Brain_Cerebellum | rs6875844 | G | -0.7140573263168335 | risk allele G decreases intron usage of chr5:151452216-151458788(+) |
| lbd | Brain_Cerebellum | rs6875844 | G | 1.1780455112457275 | risk allele G increases intron usage of chr5:151452491-151458788(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBM41, G3BP2, HNRNPLL, MATR3, FUS, CELF5, PABPC3, IGF2BP2, YTHDC1, PTBP2
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), AGO1(2), AGO2(2), CELF4(2), CELF5(2), CELF6(2), CPEB1(2), CPEB4(2), CSTF2(2), DAZAP1(2), DDX19B(2), DDX58(2), EIF4A3(2), EIF4B(2), ELAVL1(2), ELAVL3(2), ELAVL4(2), ENOX1(2), ESRP1(2)

## 5. Interpretation
SLC36A1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
