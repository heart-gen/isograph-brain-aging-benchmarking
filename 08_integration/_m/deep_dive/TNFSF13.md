# TNFSF13 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.02 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF n/a  ·  missense o/e n/a

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Caudate_basal_ganglia | 0.02 | no | — | no | — | risk allele G decreases TNFSF13 expression (gene-level; no intron) |
| als | sQTL | Brain_Putamen_basal_ganglia | 0.02 | no | no | no | — | risk allele G increases usage of junction chr17:7559297-7559846(+) (ENST00000380535.8; not |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Caudate_basal_ganglia | rs3803800 | G | -0.218720480799675 | risk allele G decreases expression of TNFSF13 |
| als | Brain_Putamen_basal_ganglia | rs3803800 | G | 0.6242228150367737 | risk allele G increases intron usage of chr17:7559297-7559846(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMY1A1, NOVA2, PABPC1, PUM1, TRA2A
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), AGO2(1), CELF4(1), CELF5(1), CMTR1(1), CPEB1(1), CSTF2(1), DAZAP1(1), G3BP2(1), IGF2BP1(1), IGHMBP2(1), KHDRBS2(1), KHDRBS3(1), LIN28A(1), NELFE(1), NONO(1), NOVA2(1), PABPC1(1), PABPC3(1), PABPC4(1)

## 5. Interpretation
TNFSF13 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
