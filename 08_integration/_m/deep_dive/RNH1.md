# RNH1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.116  ·  missense o/e 1.05

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.02 | no | no | no | — | risk allele A decreases usage of junction chr11:502181-504427(-) (ENST00000529368.5; not i |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.02 | no | no | no | — | risk allele A decreases usage of junction chr11:502181-506609(-) (ENST00000438658.6; not i |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.02 | no | no | no | — | risk allele A decreases usage of junction chr11:502181-507113(-) (ENST00000397604.7,ENST00 |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.02 | no | no | no | — | risk allele A increases usage of junction chr11:502249-504824(-) (ENST00000354420.7,ENST00 |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.02 | no | no | no | — | risk allele A increases usage of junction chr11:504996-506609(-) (ENST00000397615.6; not i |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.02 | no | no | no | — | risk allele A increases usage of junction chr11:504996-507113(-) (ENST00000354420.7,ENST00 |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Putamen_basal_ganglia | rs67092853 | A | -0.8540651202201843 | risk allele A decreases intron usage of chr11:502181-504427(-) |
| scz | Brain_Putamen_basal_ganglia | rs67092853 | A | -1.1529988050460815 | risk allele A decreases intron usage of chr11:502181-506609(-) |
| scz | Brain_Putamen_basal_ganglia | rs67092853 | A | -1.3136303424835205 | risk allele A decreases intron usage of chr11:502181-507113(-) |
| scz | Brain_Putamen_basal_ganglia | rs67092853 | A | 1.1085492372512817 | risk allele A increases intron usage of chr11:502249-504824(-) |
| scz | Brain_Putamen_basal_ganglia | rs67092853 | A | 1.1551642417907715 | risk allele A increases intron usage of chr11:504996-506609(-) |
| scz | Brain_Putamen_basal_ganglia | rs67092853 | A | 0.7876113057136536 | risk allele A increases intron usage of chr11:504996-507113(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZRANB2, NONO
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ACO1(1), AGO2(1), CELF4(1), CELF5(1), CELF6(1), CNOT4(1), CPEB1(1), CPEB4(1), CSTF2(1), ESRP1(1), FXR1(1), FXR2(1), GRSF1(1), HNRNPAB(1), HNRNPK(1), HNRNPLL(1), HNRNPM(1), HNRNPU(1), IGHMBP2(1)

## 5. Interpretation
RNH1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
