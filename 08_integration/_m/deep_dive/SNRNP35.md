# SNRNP35 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.971  ·  missense o/e 0.87

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr12:123458216-123458964(+) (ENST00000529904.2; |
| pd | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T increases usage of junction chr12:123458216-123465538(+) (ENST00000526639.3; |
| pd | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T increases usage of junction chr12:123458216-123471092(+) (ENST00000527158.2; |
| pd | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr12:123459900-123465538(+) (ENST00000412157.2; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Spinal_cord_cervical_c-1 | rs4479070 | T | -0.5112804770469666 | risk allele T decreases intron usage of chr12:123458216-123458964(+) |
| pd | Brain_Spinal_cord_cervical_c-1 | rs4479070 | T | 0.4988984167575836 | risk allele T increases intron usage of chr12:123458216-123465538(+) |
| pd | Brain_Spinal_cord_cervical_c-1 | rs4479070 | T | 0.4425353407859802 | risk allele T increases intron usage of chr12:123458216-123471092(+) |
| pd | Brain_Spinal_cord_cervical_c-1 | rs4479070 | T | -0.6075859665870667 | risk allele T decreases intron usage of chr12:123459900-123465538(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPC, HNRNPA0, KHDRBS1, DAZAP1, HNRNPU, HNRNPK, SYNCRIP, PABPC1, NUDT21, QKI, U2AF2, ELAVL3
- **All switched-motif RBPs (recurrence across regions):** AGO1(2), AGO2(2), AKAP1(2), CELF4(2), CELF5(2), CPEB1(2), CPEB2(2), CPEB4(2), CSTF2(2), DAZAP1(2), DDX19B(2), EIF4A3(2), EIF4B(2), ELAVL3(2), ENOX1(2), ESRP1(2), FXR1(2), G3BP1(2), GRSF1(2), HNRNPA0(2)

## 5. Interpretation
SNRNP35 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
