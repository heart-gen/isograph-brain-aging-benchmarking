# TTC19 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.049  ·  missense o/e 1.06

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Putamen_basal_ganglia | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr17:16000245-16001915(+) (ENST00000261647.10,E |
| pd | sQTL | Brain_Putamen_basal_ganglia | 0.01 | no | no | no | — | risk allele T increases usage of junction chr17:16000382-16001915(+) (ENST00000466729.5,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Putamen_basal_ganglia | rs11656046 | T | -0.49589890241622925 | risk allele T decreases intron usage of chr17:16000245-16001915(+) |
| pd | Brain_Putamen_basal_ganglia | rs11656046 | T | 0.47326529026031494 | risk allele T increases intron usage of chr17:16000382-16001915(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, RBM6, SNRPB2, PPRC1, G3BP2, ENOX1, ZC3H10, IGF2BP1, CNOT4, HNRNPA3, CELF6, PABPC3, PABPC5
- **All switched-motif RBPs (recurrence across regions):** A1CF(7), ACO1(7), AGO2(7), AKAP1(7), CELF6(7), CNOT4(7), CPEB2(7), CPEB4(7), DHX58(7), CSTF2(7), DDX19B(7), DDX58(7), EIF4B(7), EIF4A3(7), ENOX1(7), ESRP1(7), ZNF638(7), ESRP2(7), FXR2(7), G3BP1(7)

## 5. Interpretation
TTC19 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
