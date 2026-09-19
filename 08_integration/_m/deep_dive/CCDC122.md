# CCDC122 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.03 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.535  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Hippocampus | 0.01 | no | — | yes | — | risk allele G decreases CCDC122 expression (gene-level; no intron) |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.03 | no | — | yes | — | risk allele G increases usage of junction chr13:43874927-43879418(-) (unmapped transcript; |
| scz | sQTL | Brain_Putamen_basal_ganglia | 0.03 | no | no | yes | — | risk allele G decreases usage of junction chr13:43874927-43879631(-) (ENST00000444614.8,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Hippocampus | rs12873099 | G | -0.28592362999916077 | risk allele G decreases expression of CCDC122 |
| scz | Brain_Putamen_basal_ganglia | rs12428350 | G | 0.582214891910553 | risk allele G increases intron usage of chr13:43874927-43879418(-) |
| scz | Brain_Putamen_basal_ganglia | rs12428350 | G | -0.5633746385574341 | risk allele G decreases intron usage of chr13:43874927-43879631(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ADAR(1), AGO2(1), CELF4(1), CELF5(1), CELF6(1), CMTR1(1), CNOT4(1), CPEB1(1), CPEB2(1), CPEB4(1), CSTF2(1), DAZAP1(1), DDX19B(1), DDX58(1), DHX9(1), EIF4A3(1), ELAVL3(1), ERI1(1), ESRP2(1)

## 5. Interpretation
CCDC122 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
