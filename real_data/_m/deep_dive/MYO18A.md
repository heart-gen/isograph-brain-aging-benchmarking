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
- **All switched-motif RBPs (recurrence across regions):** AKAP1(5), CPEB2(5), DDX19B(5), ESRP1(5), ENOX1(5), EIF4A3(5), HNRNPCL1(5), G3BP2(5), G3BP1(5), PABPC4(5), MATR3(5), KHDRBS2(5), IGF2BP2(5), IGF2BP1(5), ZC3H10(5), RBMS1(5), RBMS3(5), SNRPB2(5), PABPC5(5), PABPC3(5)

## 5. Interpretation
MYO18A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
