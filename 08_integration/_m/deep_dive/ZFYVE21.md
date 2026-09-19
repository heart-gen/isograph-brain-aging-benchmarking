# ZFYVE21 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.965  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.04 | no | no | no | — | risk allele G decreases usage of junction chr14:103715979-103728908(+) (ENST00000556610.5; |
| scz | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.04 | no | no | no | — | risk allele G decreases usage of junction chr14:103729182-103732620(+) (ENST00000311141.7; |
| scz | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.04 | no | no | no | — | risk allele G increases usage of junction chr14:103729846-103732620(+) (ENST00000216602.10 |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs6576007 | G | -0.30198490619659424 | risk allele G decreases intron usage of chr14:103715979-103728908(+) |
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs6576007 | G | -0.38509470224380493 | risk allele G decreases intron usage of chr14:103729182-103732620(+) |
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs6576007 | G | 0.3866266906261444 | risk allele G increases intron usage of chr14:103729846-103732620(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, HNRNPA0, HNRNPC, ELAVL3, KHDRBS1, HNRNPD, SRSF5, FMR1, CPEB1, PCBP1, ADAR, PABPC1, CPEB4, U2AF2, PPIE
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), ADAR(5), AGO2(5), AKAP1(5), CELF4(5), CELF5(5), CELF6(5), CNOT4(5), CPEB1(5), CPEB2(5), CSTF2(5), DDX19B(5), EIF4A3(5), DHX58(5), G3BP1(5), ESRP1(5), EIF4B(5), HNRNPCL1(5), HNRNPAB(5)

## 5. Interpretation
ZFYVE21 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
