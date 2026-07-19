# PPIL2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.03 (Brain_Hypothalamus)
- **Constraint:** LOEUF 0.791  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele T increases PPIL2 expression (gene-level; no intron) |
| scz | sQTL | Brain_Hypothalamus | 0.03 | no | no | yes | — | risk allele T increases usage of junction chr22:21695496-21696031(+) (ENST00000406385.1,EN |
| scz | sQTL | Brain_Hypothalamus | 0.03 | no | no | yes | — | risk allele T decreases usage of junction chr22:21695496-21696747(+) (ENST00000335025.12,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs12484060 | T | 0.3390215337276459 | risk allele T increases expression of PPIL2 |
| scz | Brain_Hypothalamus | rs5999098 | T | 0.48395127058029175 | risk allele T increases intron usage of chr22:21695496-21696031(+) |
| scz | Brain_Hypothalamus | rs5999098 | T | -0.8351643681526184 | risk allele T decreases intron usage of chr22:21695496-21696747(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** MATR3, PABPC5, YBX2, ESRP1, ZFP36L2, IGF2BP2, SNRNP70, PABPC3, RBM46, PABPC4, RBMS1, SART3, CELF5, CNOT4
- **All switched-motif RBPs (recurrence across regions):** AKAP1(4), CELF5(4), CNOT4(4), DHX58(4), ESRP1(4), IGF2BP2(4), FXR2(4), MATR3(4), PABPC3(4), YBX2(4), PABPC5(4), PHAX(4), PTBP2(4), RBM46(4), RBM25(4), RBMS1(4), SAMD4A(4), ZC3H10(4), SART3(4), SNRNP70(4)

## 5. Interpretation
PPIL2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
