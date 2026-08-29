# CDIP1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.132  ·  missense o/e 1.00

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.03 | yes | yes | yes | no annotated structural change | risk allele C decreases usage of junction chr16:4514664-4538325(-) (ENST00000563332.6,ENST |
| scz | sQTL | Brain_Caudate_basal_ganglia | 0.03 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr16:4514664-4538702(-) (ENST00000562334.5,ENST |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs4786493 | C | -0.4410758316516876 | risk allele C decreases intron usage of chr16:4514664-4538325(-) |
| scz | Brain_Caudate_basal_ganglia | rs4786493 | C | 0.43002358078956604 | risk allele C increases intron usage of chr16:4514664-4538702(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** KHDRBS1, CPEB1, CPEB4, SNRPA, NUDT21, U2AF2, HNRNPK
- **All switched-motif RBPs (recurrence across regions):** AGO2(8), CNOT4(8), CPEB1(8), CPEB2(8), CPEB4(8), DDX19B(8), DHX58(8), EIF4A3(8), FXR2(8), HNRNPAB(8), HNRNPK(8), HNRNPLL(8), HNRNPM(8), HNRNPU(8), IGF2BP2(8), IGF2BP3(8), IGHMBP2(8), KHDRBS1(8), KHDRBS2(8), KHDRBS3(8)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for CDIP1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
CDIP1 (cell-death-inducing p53 target) colocalizes as a two-event splicing-led switch in schizophrenia; its isoform biology in SCZ is uncharacterized -- a novel candidate.
