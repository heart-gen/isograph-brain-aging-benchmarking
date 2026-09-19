# MYO18A — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.397  ·  missense o/e 0.81

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | — | no | — | risk allele G increases usage of junction chr17:29074914-29077028(-) (unmapped transcript; |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr17:29074914-29082316(-) (ENST00000527312.5,EN |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele G increases usage of junction chr17:29081000-29082316(-) (ENST00000636100.1; n |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr17:29082438-29085604(-) (ENST00000527372.7; n |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele G increases usage of junction chr17:29082438-29086438(-) (ENST00000527312.5,EN |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr17:29085648-29086438(-) (ENST00000527372.7; n |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele G increases usage of junction chr17:29086577-29086936(-) (ENST00000527372.7,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Caudate_basal_ganglia | rs4965973 | G | 0.5108667016029358 | risk allele G increases intron usage of chr17:29074914-29077028(-) |
| als | Brain_Caudate_basal_ganglia | rs4965973 | G | -0.5754358768463135 | risk allele G decreases intron usage of chr17:29074914-29082316(-) |
| als | Brain_Caudate_basal_ganglia | rs4965973 | G | 0.39853960275650024 | risk allele G increases intron usage of chr17:29081000-29082316(-) |
| als | Brain_Caudate_basal_ganglia | rs4965973 | G | -0.3535462021827698 | risk allele G decreases intron usage of chr17:29082438-29085604(-) |
| als | Brain_Caudate_basal_ganglia | rs4965973 | G | 0.25569310784339905 | risk allele G increases intron usage of chr17:29082438-29086438(-) |
| als | Brain_Caudate_basal_ganglia | rs4965973 | G | -0.36540544033050537 | risk allele G decreases intron usage of chr17:29085648-29086438(-) |
| als | Brain_Caudate_basal_ganglia | rs4965973 | G | 0.4223220646381378 | risk allele G increases intron usage of chr17:29086577-29086936(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, G3BP1, RBM6, RBMS3, G3BP2, ENOX1, IGF2BP1, CNOT4, ZC3H10, HNRNPA3, CELF4, CELF5, CELF6, MATR3, PABPC3
- **All switched-motif RBPs (recurrence across regions):** AKAP1(9), CPEB4(9), DDX19B(9), CELF6(9), CNOT4(9), EIF4A3(9), IGHMBP2(9), KHDRBS2(9), MATR3(9), MSI1(9), PABPC3(9), PABPC5(9), HNRNPU(9), HNRNPCL1(9), IGF2BP1(9), G3BP2(9), G3BP1(9), FXR2(9), ESRP2(9), ESRP1(9)

## 5. Interpretation
MYO18A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
