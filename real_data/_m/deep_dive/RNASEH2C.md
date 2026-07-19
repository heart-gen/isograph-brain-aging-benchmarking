# RNASEH2C — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.05 (Brain_Anterior_cingulate_cortex_BA24)
- **Constraint:** LOEUF 1.688  ·  missense o/e 1.28

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Hypothalamus | 0.02 | no | — | yes | — | risk allele A decreases RNASEH2C expression (gene-level; no intron) |
| scz | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.05 | no | no | yes | — | risk allele T increases usage of junction chr11:65715317-65719681(-) (ENST00000644198.1; n |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Hypothalamus | rs7944860 | A | -0.328283429145813 | risk allele A decreases expression of RNASEH2C |
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs4930156 | T | 0.5024814605712891 | risk allele T increases intron usage of chr11:65715317-65719681(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** CPEB1(3), CPEB4(3), DHX58(3), ENOX1(3), ESRP2(3), FXR2(3), HNRNPAB(3), HNRNPA3(3), SUPV3L1(3), PUM2(3), U2AF2(3), RBM24(3), SNRPB2(3), DDX58(2), G3BP1(2), KHDRBS1(2), PHAX(1), HNRNPLL(1)

## 5. Interpretation
RNASEH2C has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
