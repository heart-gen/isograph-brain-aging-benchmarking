# TMEM63A — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.839  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | yes | — | risk allele C decreases usage of junction chr1:225865967-225866574(-) (ENST00000366835.8;  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Spinal_cord_cervical_c-1 | rs168617 | C | -0.5349173545837402 | risk allele C decreases intron usage of chr1:225865967-225866574(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, PABPN1, HNRNPC, HNRNPA0, PUM2, CELF6, RBM25, KHDRBS1, SUPV3L1, IGF2BP2, DAZAP1, HNRNPU, PIWIL1, HNRNPK, SYNCRIP
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), ADAR(3), AGO1(3), AGO2(3), CELF4(3), AKAP1(3), CELF5(3), CELF6(3), CSTF2(3), CNOT4(3), CPEB1(3), CPEB2(3), CPEB4(3), EIF4A3(3), DHX58(3), DDX19B(3), HNRNPA0(3), FXR1(3), G3BP1(3), EIF4B(3)

## 5. Interpretation
TMEM63A has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
