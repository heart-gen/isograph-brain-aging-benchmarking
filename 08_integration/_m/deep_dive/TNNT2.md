# TNNT2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Amygdala)
- **Constraint:** LOEUF n/a  ·  missense o/e n/a

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Amygdala | 0.01 | no | no | no | — | risk allele C increases usage of junction chr1:201367806-201368162(-) (ENST00000367318.10, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Amygdala | rs16848211 | C | 2.403996229171753 | risk allele C increases intron usage of chr1:201367806-201368162(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** —

## 5. Interpretation
TNNT2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
