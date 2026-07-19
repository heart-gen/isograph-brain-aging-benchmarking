# GABBR2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.381  ·  missense o/e 0.39

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | no | — | risk allele G decreases usage of junction chr9:98299353-98303241(-) (ENST00000259455.4,ENS |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | no | — | risk allele G increases usage of junction chr9:98302896-98303241(-) (ENST00000636575.1; no |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs7869257 | G | -0.6419111490249634 | risk allele G decreases intron usage of chr9:98299353-98303241(-) |
| scz | Brain_Cerebellum | rs7869257 | G | 0.5653102993965149 | risk allele G increases intron usage of chr9:98302896-98303241(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PPIE, ELAVL4, KHDRBS1, HNRNPD, SYNCRIP, RNASEL, ELAVL3, OAS1, IGF2BP3, TIAL1, U2AF2, ELAVL2, PABPC1, KHDRBS3, RBMX
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO1(4), AGO2(4), CPEB1(4), CPEB4(4), DHX58(4), DDX58(4), DDX19B(4), CSTF2(4), HNRNPLL(4), HNRNPD(4), HNRNPCL1(4), HNRNPA3(4), EIF4A3(4), DHX9(4), ELAVL4(4), ELAVL3(4), ESRP2(4)

## 5. Interpretation
GABBR2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
