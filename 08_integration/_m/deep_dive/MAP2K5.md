# MAP2K5 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.201  ·  missense o/e 0.87

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | — | no | — | risk allele G decreases usage of junction chr15:67563350-67585890(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs12905509 | G | -0.32863107323646545 | risk allele G decreases intron usage of chr15:67563350-67585890(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS3, SNRPB2, PPRC1, RBMS1, FUS, PABPC4, G3BP1, PABPC3, DHX9, YTHDC1, PTBP2, ADAR
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), ADAR(3), AGO1(3), AGO2(3), ANKHD1(3), CPEB1(3), CSTF2(3), DDX19B(3), DDX58(3), DHX9(3), DHX58(3), EIF4A3(3), EIF4B(3), ZRANB2(3), ELAVL3(3), ESRP1(3), G3BP1(3), IFIH1(3), HNRNPU(3), HNRNPCL1(3)

## 5. Interpretation
MAP2K5 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
