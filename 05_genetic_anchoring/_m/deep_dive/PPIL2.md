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
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, RBMS1, G3BP1, ZC3H10, FXR2, ESRP1, YTHDC1, HNRNPA3, CNOT4, CPEB2, ENOX1, PABPC3, PABPN1, PABPC5, RBM24
- **All switched-motif RBPs (recurrence across regions):** ESRP1(6), ENOX1(6), CPEB2(6), CNOT4(6), HNRNPA3(6), FXR2(6), PHAX(6), PABPN1(6), PTBP2(6), PUM2(6), IGF2BP2(6), HNRNPCL1(6), ZC3H10(6), YTHDC1(6), SAMD4A(6), RBMS1(6), RBM41(6), RBM46(6), PABPC5(6), PABPC3(6)

## 5. Interpretation
PPIL2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
