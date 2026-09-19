# RAD51C — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** AD,SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.06 (Brain_Cortex)
- **Constraint:** LOEUF 1.134  ·  missense o/e 0.89

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cortex | 0.06 | no | no | no | — | risk allele A decreases usage of junction chr17:58692959-58694931(+) (ENST00000697692.1; n |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | no | — | direction unresolved |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cortex | rs2526377 | A | -0.49035006761550903 | risk allele A decreases intron usage of chr17:58692959-58694931(+) |
| scz | Brain_Cerebellum | rs11650105 | C | nan | unresolved (variant not in GTEx signif_pairs or allele mismatch) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SNRPB2, PPRC1, HNRNPLL, RBMS1, IGF2BP2, YTHDC1, PTBP2
- **All switched-motif RBPs (recurrence across regions):** AGO2(1), AKAP1(1), CMTR1(1), CPEB2(1), CPEB4(1), CSTF2(1), DDX19B(1), DDX58(1), EIF4A3(1), ELAVL3(1), ESRP1(1), FMR1(1), FXR1(1), HNRNPAB(1), HNRNPCL1(1), HNRNPK(1), HNRNPLL(1), HNRNPM(1), HNRNPU(1), IGF2BP2(1)

## 5. Interpretation
RAD51C has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
