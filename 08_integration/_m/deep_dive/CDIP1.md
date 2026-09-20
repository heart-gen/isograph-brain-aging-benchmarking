# CDIP1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible · replicates in BrainSeq

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.132  ·  missense o/e 1.00

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.03 | no | no | yes | — | risk allele C decreases usage of junction chr16:4514664-4538325(-) (ENST00000563332.6,ENST |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.03 | no | no | yes | — | risk allele C increases usage of junction chr16:4514664-4538702(-) (ENST00000562334.5,ENST |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs4786493 | C | -0.4410758316516876 | risk allele C decreases intron usage of chr16:4514664-4538325(-) |
| scz | Brain_Caudate_basal_ganglia | rs4786493 | C | 0.43002358078956604 | risk allele C increases intron usage of chr16:4514664-4538702(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS3, RBMY1A1, IGF2BP3, HNRNPLL, RBM14, IGF2BP2, PTBP2, RBM46
- **All switched-motif RBPs (recurrence across regions):** AGO2(6), CNOT4(6), CPEB1(6), CPEB2(6), CPEB4(6), DDX19B(6), DHX58(6), EIF4A3(6), FXR2(6), HNRNPAB(6), HNRNPCL1(6), HNRNPLL(6), HNRNPM(6), IGF2BP2(6), IGF2BP3(6), IGHMBP2(6), KHDRBS1(6), KHDRBS2(6), KHDRBS3(6), LIN28A(6)

## 5. Interpretation
CDIP1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.

## 6. Literature (known isoform biology)
CDIP1 (cell-death-inducing p53 target) was nominated as a two-event splicing-led switch in schizophrenia; its isoform biology is uncharacterised. Not in the current anchored set.

_Curation: novel_candidate._ No established disease-specific isoform biology was found. That is a statement about the literature, not about the evidence here: an under-characterised switch is what this method is built to surface, and this row is a nomination rather than a null result.
