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
- **Switched *and* module-enriched (q<0.05) RBPs:** ELAVL4, KHDRBS1, TIAL1, RNASEL, PPIE, ZFP36, OAS1, HNRNPC, TIA1, IGF2BP3, HNRNPD, ELAVL2, PABPC1, RBMX, ELAVL3
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), ADAR(6), AGO1(6), AGO2(6), ANKHD1(6), CPEB2(6), CPEB1(6), CNOT4(6), CPEB4(6), DDX58(6), DDX19B(6), CSTF2(6), EIF4A3(6), DHX9(6), DHX58(6), FUS(6), ESRP1(6), ESRP2(6), ENOX1(6)

## 5. Interpretation
GABBR2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
