# DOC2A — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** AD,SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.05 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.732  ·  missense o/e 0.84

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Cerebellar_Hemisphere | 0.04 | no | — | no | — | risk allele C decreases DOC2A expression (gene-level; no intron) |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.05 | no | no | no | — | risk allele C decreases usage of junction chr16:30007090-30007173(-) (ENST00000350119.9,EN |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.05 | no | — | no | — | risk allele C increases usage of junction chr16:30007090-30007179(-) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs3814883 | C | -0.7284601330757141 | risk allele C decreases intron usage of chr16:30007090-30007173(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs3814883 | C | 0.46809905767440796 | risk allele C increases intron usage of chr16:30007090-30007179(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, ELAVL3, PPIE, CPEB4, HNRNPDL, HNRNPAB, LIN28A, HNRNPLL, DHX58, RBFOX2, TIA1, RBM4, SNRPA, ELAVL2, U2AF2
- **All switched-motif RBPs (recurrence across regions):** ANKHD1(5), RALY(5), CELF4(5), CELF5(5), CSTF2(5), CPEB4(5), DDX19B(5), ESRP2(5), HNRNPLL(5), IGF2BP2(5), KHDRBS2(5), HNRNPCL1(5), NELFE(5), PTBP2(5), PIWIL1(5), PABPC5(5), PABPC3(5), RBM4(5), RBM46(5), RBM6(5)

## 5. Interpretation
DOC2A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
