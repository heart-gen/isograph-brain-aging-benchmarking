# NMRAL1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF n/a  ·  missense o/e n/a

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.03 | no | no | no | — | risk allele A increases usage of junction chr16:4466402-4469504(-) (ENST00000571291.5; not |
| scz | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.03 | no | no | no | — | risk allele A increases usage of junction chr16:4471449-4474093(-) (ENST00000573520.5,ENST |
| scz | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.03 | no | no | no | — | risk allele A increases usage of junction chr16:4474166-4474423(-) (ENST00000574425.5,ENST |
| scz | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.03 | no | no | no | — | risk allele A decreases usage of junction chr16:4474166-4474554(-) (ENST00000283429.11,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Spinal_cord_cervical_c-1 | rs4786510 | A | 0.5825082063674927 | risk allele A increases intron usage of chr16:4466402-4469504(-) |
| scz | Brain_Spinal_cord_cervical_c-1 | rs4786510 | A | 0.5299429297447205 | risk allele A increases intron usage of chr16:4471449-4474093(-) |
| scz | Brain_Spinal_cord_cervical_c-1 | rs4786510 | A | 1.073625922203064 | risk allele A increases intron usage of chr16:4474166-4474423(-) |
| scz | Brain_Spinal_cord_cervical_c-1 | rs4786510 | A | -1.1471813917160034 | risk allele A decreases intron usage of chr16:4474166-4474554(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, SAMD4A, PABPC5
- **All switched-motif RBPs (recurrence across regions):** AGO2(1), CELF1(1), CELF2(1), CPEB2(1), DHX58(1), ELAVL3(1), ESRP2(1), FXR2(1), G3BP1(1), HNRNPA0(1), HNRNPA3(1), HNRNPAB(1), HNRNPC(1), HNRNPD(1), HNRNPK(1), HNRNPLL(1), HNRNPU(1), IGHMBP2(1), KHDRBS1(1), KHDRBS2(1)

## 5. Interpretation
NMRAL1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
