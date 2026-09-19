# LRRC63 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 1.303  ·  missense o/e 0.93

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | — | no | — | risk allele A decreases usage of junction chr13:46212259-46228665(+) (unmapped transcript; |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele A increases usage of junction chr13:46246625-46257136(+) (ENST00000676025.1; n |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | — | no | — | risk allele A decreases usage of junction chr13:46246625-46266733(+) (unmapped transcript; |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | — | no | — | risk allele A decreases usage of junction chr13:46246625-46276590(+) (unmapped transcript; |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele A increases usage of junction chr13:46254119-46257136(+) (ENST00000674570.1,EN |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele A increases usage of junction chr13:46254489-46257136(+) (ENST00000674850.1,EN |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | — | no | — | risk allele A increases usage of junction chr13:46257273-46266733(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Spinal_cord_cervical_c-1 | rs1749644 | A | -0.6625434756278992 | risk allele A decreases intron usage of chr13:46212259-46228665(+) |
| als | Brain_Spinal_cord_cervical_c-1 | rs1749644 | A | 0.8294452428817749 | risk allele A increases intron usage of chr13:46246625-46257136(+) |
| als | Brain_Spinal_cord_cervical_c-1 | rs1749644 | A | -1.0238467454910278 | risk allele A decreases intron usage of chr13:46246625-46266733(+) |
| als | Brain_Spinal_cord_cervical_c-1 | rs1749644 | A | -0.5055336356163025 | risk allele A decreases intron usage of chr13:46246625-46276590(+) |
| als | Brain_Spinal_cord_cervical_c-1 | rs1749644 | A | 0.9940190315246582 | risk allele A increases intron usage of chr13:46254119-46257136(+) |
| als | Brain_Spinal_cord_cervical_c-1 | rs1749644 | A | 0.7252926230430603 | risk allele A increases intron usage of chr13:46254489-46257136(+) |
| als | Brain_Spinal_cord_cervical_c-1 | rs1749644 | A | 0.5919929146766663 | risk allele A increases intron usage of chr13:46257273-46266733(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SRSF10
- **All switched-motif RBPs (recurrence across regions):** AGO2(1), AKAP1(1), CELF4(1), CELF5(1), CELF6(1), CSTF2(1), ENOX1(1), G3BP2(1), HNRNPAB(1), HNRNPLL(1), KHDRBS2(1), KHDRBS3(1), MSI1(1), NELFE(1), PABPN1(1), PUM2(1), RBM24(1), RBM3(1), RBM6(1), RBMS1(1)

## 5. Interpretation
LRRC63 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
