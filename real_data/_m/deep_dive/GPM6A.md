# GPM6A — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.21 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.362  ·  missense o/e 0.61

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.07 | no | — | yes | — | risk allele G decreases GPM6A expression (gene-level; no intron) |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.21 | no | no | yes | — | risk allele G decreases usage of junction chr4:175701767-175812191(-) (ENST00000280187.11, |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.21 | no | no | yes | — | risk allele G increases usage of junction chr4:175701767-176002309(-) (ENST00000506894.5;  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs78640136 | G | -0.27255088090896606 | risk allele G decreases expression of GPM6A |
| scz | Brain_Caudate_basal_ganglia | rs78640136 | G | -0.717190682888031 | risk allele G decreases intron usage of chr4:175701767-175812191(-) |
| scz | Brain_Caudate_basal_ganglia | rs78640136 | G | 0.7075767517089844 | risk allele G increases intron usage of chr4:175701767-176002309(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZCRB1, CPEB2, DHX9, RBMS3, ZFP36L2, DDX58, AKAP1, IFIH1, RBM42, RALY, ADAR, DHX58, PABPC5, PABPC4, ZNF346
- **All switched-motif RBPs (recurrence across regions):** ADAR(5), AGO1(5), AKAP1(5), CELF4(5), CELF5(5), CMTR1(5), CPEB1(5), CPEB2(5), CSTF2(5), DDX19B(5), DDX58(5), DHX58(5), DHX9(5), EIF4A3(5), ERI1(5), ESRP1(5), ESRP2(5), G3BP1(5), HNRNPA0(5), HNRNPAB(5)

## 5. Interpretation
GPM6A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
