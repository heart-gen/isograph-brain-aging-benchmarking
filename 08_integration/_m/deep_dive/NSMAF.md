# NSMAF — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.767  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | — | no | — | risk allele T decreases NSMAF expression (gene-level; no intron) |
| als | sQTL | Brain_Cerebellum | 0.04 | no | no | no | — | risk allele T increases usage of junction chr8:58635546-58642984(-) (ENST00000038176.8,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs75712986 | T | -0.6134388446807861 | risk allele T decreases expression of NSMAF |
| als | Brain_Cerebellum | rs75712986 | T | 1.1249518394470215 | risk allele T increases intron usage of chr8:58635546-58642984(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBM41, SNRPB2, PPRC1, HNRNPLL, RBM14, RBMS1, MATR3, CELF5, DHX9, IGF2BP2, YTHDC1, ADAR
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), ADAR(2), AGO1(2), AGO2(2), ANKHD1(2), CELF4(2), CELF5(2), CNOT4(2), CPEB2(2), CPEB4(2), CSTF2(2), DAZAP1(2), DDX19B(2), DDX58(2), DHX58(2), DHX9(2), EIF4A3(2), ELAVL3(2), ENOX1(2)

## 5. Interpretation
NSMAF has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
