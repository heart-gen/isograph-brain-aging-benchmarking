# BCKDK — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.07 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.439  ·  missense o/e 0.95

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.07 | no | — | no | — | risk allele C decreases usage of junction chr16:31115869-31117463(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Spinal_cord_cervical_c-1 | rs889555 | C | -0.5409450531005859 | risk allele C decreases intron usage of chr16:31115869-31117463(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, G3BP2, AKAP1, SNRPB2, SNRNP70
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), AGO1(2), AGO2(2), AKAP1(2), CELF1(2), CELF2(2), CELF6(2), CPEB1(2), CSTF2(2), DAZAP1(2), DHX58(2), ELAVL1(2), ELAVL3(2), ELAVL4(2), ESRP2(2), F2(2), FXR2(2), G3BP2(2), HNRNPA2B1(2)

## 5. Interpretation
BCKDK has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
