# FAM221A — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible · replicates in BrainSeq

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.554  ·  missense o/e 1.10

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C increases usage of junction chr7:23691596-23692183(+) (ENST00000429719.5,ENS |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C decreases usage of junction chr7:23691596-23698192(+) (ENST00000344962.9,ENS |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C increases usage of junction chr7:23691596-23700786(+) (ENST00000409192.7,ENS |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C decreases usage of junction chr7:23698299-23700786(+) (ENST00000344962.9,ENS |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | — | yes | — | risk allele C decreases usage of junction chr7:23698299-23700789(+) (unmapped transcript;  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs5029444 | C | 0.4050355553627014 | risk allele C increases intron usage of chr7:23691596-23692183(+) |
| scz | Brain_Caudate_basal_ganglia | rs5029444 | C | -0.8854654431343079 | risk allele C decreases intron usage of chr7:23691596-23698192(+) |
| scz | Brain_Caudate_basal_ganglia | rs5029444 | C | 1.099471092224121 | risk allele C increases intron usage of chr7:23691596-23700786(+) |
| scz | Brain_Caudate_basal_ganglia | rs5029444 | C | -1.0529088973999023 | risk allele C decreases intron usage of chr7:23698299-23700786(+) |
| scz | Brain_Caudate_basal_ganglia | rs5029444 | C | -0.5256322026252747 | risk allele C decreases intron usage of chr7:23698299-23700789(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB2, A1CF, RBMS3, ZFP36L2, RBMS1, AKAP1, RBM41, RBM42, RALY, DHX58, PABPC5, PABPC4, SART3, IGF2BP2
- **All switched-motif RBPs (recurrence across regions):** AGO2(4), AGO1(4), CELF4(4), AKAP1(4), CSTF2(4), EIF4A3(4), CELF5(4), CPEB2(4), ESRP1(4), HNRNPA2B1(4), HNRNPA0(4), G3BP1(4), ZFP36L2(4), YTHDC1(4), SRSF10(4), HNRNPAB(4), HNRNPD(4), HNRNPCL1(4), HNRNPLL(4), IFIH1(4)

## 5. Interpretation
FAM221A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
