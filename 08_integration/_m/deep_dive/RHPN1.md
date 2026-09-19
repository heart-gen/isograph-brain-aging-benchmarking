# RHPN1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.605  ·  missense o/e 1.19

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.04 | no | no | yes | — | risk allele C increases usage of junction chr8:143377455-143378269(+) (ENST00000289013.11, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs117605463 | C | 3.0571868419647217 | risk allele C increases intron usage of chr8:143377455-143378269(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), AGO1(2), AGO2(2), ELAVL1(2), ELAVL3(2), ELAVL4(2), FXR1(2), G3BP1(2), HNRNPA0(2), HNRNPD(2), HNRNPDL(2), HNRNPLL(2), HNRNPU(2), IGHMBP2(2), KHDRBS1(2), KHDRBS2(2), KHDRBS3(2), NELFE(2), NOVA2(2), NUDT21(2)

## 5. Interpretation
RHPN1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
