# MYO18A — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible · replicates in BrainSeq

- **Traits:** ALS,SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.397  ·  missense o/e 0.81

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | — | yes | — | risk allele G increases usage of junction chr17:29074914-29077028(-) (unmapped transcript; |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele G decreases usage of junction chr17:29074914-29082316(-) (ENST00000527312.5,EN |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr17:29081000-29082316(-) (ENST00000636100.1; n |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele G decreases usage of junction chr17:29082438-29085604(-) (ENST00000527372.7; n |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr17:29082438-29086438(-) (ENST00000527312.5,EN |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele G decreases usage of junction chr17:29085648-29086438(-) (ENST00000527372.7; n |
| als | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr17:29086577-29086936(-) (ENST00000527372.7,EN |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | — | no | — | risk allele A decreases usage of junction chr17:29074914-29077028(-) (unmapped transcript; |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele A increases usage of junction chr17:29074914-29082316(-) (ENST00000527312.5,EN |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele A decreases usage of junction chr17:29081000-29082316(-) (ENST00000636100.1; n |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele A increases usage of junction chr17:29082438-29085604(-) (ENST00000527372.7; n |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele A decreases usage of junction chr17:29082438-29086438(-) (ENST00000527312.5,EN |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele A increases usage of junction chr17:29085648-29086438(-) (ENST00000527372.7; n |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | no | — | risk allele A decreases usage of junction chr17:29086577-29086936(-) (ENST00000527372.7,EN |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** IGF2BP1, ESRP1, PABPC3, MATR3, G3BP1, CNOT4
- **All switched-motif RBPs (recurrence across regions):** AGO2(6), AKAP1(6), CELF6(6), CNOT4(6), DDX19B(6), EIF4A3(6), ESRP1(6), ENOX1(6), G3BP1(6), FXR2(6), ZFP36L2(6), ZNF638(6), G3BP2(6), HNRNPCL1(6), IGF2BP2(6), IGF2BP1(6), MATR3(6), KHDRBS2(6), PABPC4(6), PABPC3(6)

## 5. Interpretation
MYO18A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
