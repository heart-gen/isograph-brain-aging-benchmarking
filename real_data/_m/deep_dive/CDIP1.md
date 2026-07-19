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
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB1, PPRC1, SUPV3L1, CPEB4, IGHMBP2
- **All switched-motif RBPs (recurrence across regions):** AGO2(6), AKAP1(6), CPEB1(6), CPEB2(6), CPEB4(6), DHX58(6), EIF4A3(6), FXR2(6), G3BP1(6), HNRNPLL(6), IGF2BP2(6), IGF2BP3(6), IGHMBP2(6), NUDT21(6), PABPN1(6), PPRC1(6), PTBP2(6), RALY(6), RBFOX2(6), RBM14(6)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for CDIP1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
CDIP1 (cell-death-inducing p53 target) colocalizes as a two-event splicing-led switch in schizophrenia; its isoform biology in SCZ is uncharacterized -- a novel candidate.
