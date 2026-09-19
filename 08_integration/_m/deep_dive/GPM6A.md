# GPM6A — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · replicates in BrainSeq

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.25 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.362  ·  missense o/e 0.61

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.09 | no | — | no | — | risk allele G decreases GPM6A expression (gene-level; no intron) |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.25 | no | no | no | — | risk allele G decreases usage of junction chr4:175701767-175812191(-) (ENST00000280187.11, |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.25 | no | no | no | — | risk allele G increases usage of junction chr4:175701767-176002309(-) (ENST00000506894.5;  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs78640136 | G | -0.27255088090896606 | risk allele G decreases expression of GPM6A |
| scz | Brain_Caudate_basal_ganglia | rs78640136 | G | -0.717190682888031 | risk allele G decreases intron usage of chr4:175701767-175812191(-) |
| scz | Brain_Caudate_basal_ganglia | rs78640136 | G | 0.7075767517089844 | risk allele G increases intron usage of chr4:175701767-176002309(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS3, SNRPB2, PPRC1, ERI1, PABPC4, HNRNPLL, RBM14, MATR3, CELF5, G3BP1, DHX9, IGF2BP2, YTHDC1, PTBP2, HNRNPA3
- **All switched-motif RBPs (recurrence across regions):** ADAR(3), AGO1(3), CELF4(3), CELF5(3), CSTF2(3), CMTR1(3), CPEB1(3), CPEB2(3), DDX58(3), DDX19B(3), DHX9(3), DHX58(3), ZFP36L2(3), ZNF346(3), EIF4A3(3), ELAVL3(3), ERI1(3), ESRP1(3), ESRP2(3), FMR1(3)

## 5. Interpretation
GPM6A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
