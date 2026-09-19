# LMF1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.161  ·  missense o/e 1.16

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | no | yes | — | risk allele C increases usage of junction chr16:954666-981145(-) (ENST00000545827.6,ENST00 |
| als | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | yes | — | risk allele C decreases usage of junction chr16:979110-979601(-) (unmapped transcript; not |
| als | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | yes | — | risk allele C increases usage of junction chr16:979110-980127(-) (unmapped transcript; not |
| als | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | yes | — | risk allele C increases usage of junction chr16:979110-981145(-) (unmapped transcript; not |
| als | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | yes | — | risk allele C decreases usage of junction chr16:979758-981145(-) (unmapped transcript; not |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Frontal_Cortex_BA9 | rs11866692 | C | 0.6495015621185303 | risk allele C increases intron usage of chr16:954666-981145(-) |
| als | Brain_Frontal_Cortex_BA9 | rs11866692 | C | -1.6825580596923828 | risk allele C decreases intron usage of chr16:979110-979601(-) |
| als | Brain_Frontal_Cortex_BA9 | rs11866692 | C | 1.3831530809402466 | risk allele C increases intron usage of chr16:979110-980127(-) |
| als | Brain_Frontal_Cortex_BA9 | rs11866692 | C | 1.7107981443405151 | risk allele C increases intron usage of chr16:979110-981145(-) |
| als | Brain_Frontal_Cortex_BA9 | rs11866692 | C | -0.6874407529830933 | risk allele C decreases intron usage of chr16:979758-981145(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, G3BP1, G3BP2, CNOT4, IGF2BP1, ZC3H10, HNRNPA3, FXR1, SNRNP70, PABPC3, PABPC5, YTHDC1, EIF4B, CPEB2, IGF2BP2
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), CPEB4(5), CPEB2(5), ELAVL3(5), DHX58(5), CSTF2(5), ELAVL1(5), PABPC5(5), HNRNPLL(5), HNRNPCL1(5), HNRNPAB(5), HNRNPA3(5), FXR2(5), G3BP1(5), SART3(5), SF1(5), TIAL1(5), SNRPA(5), TRA2A(5), U2AF2(5)

## 5. Interpretation
LMF1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
