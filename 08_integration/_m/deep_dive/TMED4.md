# TMED4 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.06 (Brain_Hypothalamus)
- **Constraint:** LOEUF 1.294  ·  missense o/e 1.02

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.02 | no | — | yes | — | risk allele G decreases TMED4 expression (gene-level; no intron) |
| scz | sQTL | Brain_Hypothalamus | 0.06 | no | no | yes | — | risk allele G increases usage of junction chr7:44579628-44581093(-) (ENST00000457408.7; no |
| scz | sQTL | Brain_Hypothalamus | 0.06 | no | no | yes | — | risk allele G decreases usage of junction chr7:44579628-44581449(-) (ENST00000289577.10; n |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Spinal_cord_cervical_c-1 | rs217384 | G | -0.2239401936531067 | risk allele G decreases expression of TMED4 |
| scz | Brain_Hypothalamus | rs217384 | G | 0.4946306645870209 | risk allele G increases intron usage of chr7:44579628-44581093(-) |
| scz | Brain_Hypothalamus | rs217384 | G | -0.6302598118782043 | risk allele G decreases intron usage of chr7:44579628-44581449(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPC, HNRNPA0, IGF2BP3, KHDRBS1, ELAVL4, DAZAP1, HNRNPU, HNRNPK, SYNCRIP, PABPC1, ELAVL1, ZFP36, NUDT21, RBMS1, QKI
- **All switched-motif RBPs (recurrence across regions):** AGO1(5), AGO2(5), AKAP1(5), CPEB1(5), CPEB2(5), CPEB4(5), CSTF2(5), DAZAP1(5), DDX19B(5), DHX58(5), EIF4A3(5), ELAVL1(5), ELAVL3(5), ELAVL4(5), ERI1(5), G3BP1(5), GRSF1(5), HNRNPA0(5), HNRNPA1L2(5), HNRNPA3(5)

## 5. Interpretation
TMED4 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
