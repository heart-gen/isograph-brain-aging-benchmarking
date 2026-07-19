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
- **Switched *and* module-enriched (q<0.05) RBPs:** YTHDC1, SUPV3L1, SNRPB2, DHX58, HNRNPAB, AGO2, CELF4, PABPC4, SYNCRIP, PUM2, CELF5
- **All switched-motif RBPs (recurrence across regions):** ANKHD1(3), CELF4(3), CELF6(3), CELF5(3), CPEB4(3), DDX19B(3), SRSF10(3), DHX58(3), FXR2(3), G3BP1(3), HNRNPAB(3), HNRNPA1L2(3), KHDRBS2(3), HNRNPLL(3), PABPN1(3), PIWIL1(3), SUPV3L1(3), ZRANB2(3), YTHDC1(3), SNRNP70(3)

## 5. Interpretation
DOC2A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
