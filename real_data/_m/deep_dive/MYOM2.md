# MYOM2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.464  ·  missense o/e 1.39

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele G decreases usage of junction chr8:2124217-2129127(+) (ENST00000262113.9,ENST0 |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele G increases usage of junction chr8:2127881-2129127(+) (ENST00000612167.4,ENST0 |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs12681998 | G | -0.6944310069084167 | risk allele G decreases intron usage of chr8:2124217-2129127(+) |
| scz | Brain_Caudate_basal_ganglia | rs12681998 | G | 0.6909205913543701 | risk allele G increases intron usage of chr8:2127881-2129127(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** IGF2BP2, PPRC1, SUPV3L1, SNRPB2, HNRNPAB, AGO1, ZFP36L2, RBM28, CELF4, A1CF
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), ADAR(3), AGO1(3), CELF4(3), CPEB1(3), CPEB4(3), CSTF2(3), DHX9(3), DHX58(3), ESRP2(3), HNRNPA1L2(3), HNRNPA3(3), HNRNPA2B1(3), TRA2A(3), SRSF11(3), HNRNPAB(3), HNRNPLL(3), IGF2BP1(3), IGF2BP2(3)

## 5. Interpretation
MYOM2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
