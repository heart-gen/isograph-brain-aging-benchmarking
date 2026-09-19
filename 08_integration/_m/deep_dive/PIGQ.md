# PIGQ — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.275  ·  missense o/e 1.05

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Putamen_basal_ganglia | 0.01 | yes | yes | yes | no annotated structural change | risk allele G increases usage of junction chr16:580972-582248(+) (ENST00000321878.10,ENST0 |
| als | sQTL | Brain_Putamen_basal_ganglia | 0.01 | yes | yes | yes | no annotated structural change | risk allele G decreases usage of junction chr16:580972-582883(+) (ENST00000026218.9; match |
| als | sQTL | Brain_Putamen_basal_ganglia | 0.01 | yes | yes | yes | no annotated structural change | risk allele G increases usage of junction chr16:582309-582883(+) (ENST00000321878.10,ENST0 |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Putamen_basal_ganglia | rs2891650 | G | 0.5695927143096924 | risk allele G increases intron usage of chr16:580972-582248(+) |
| als | Brain_Putamen_basal_ganglia | rs2891650 | G | -1.0913244485855103 | risk allele G decreases intron usage of chr16:580972-582883(+) |
| als | Brain_Putamen_basal_ganglia | rs2891650 | G | 0.6233781576156616 | risk allele G increases intron usage of chr16:582309-582883(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PPRC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), AGO1(2), CMTR1(2), CPEB1(2), CPEB2(2), CPEB4(2), CSTF2(2), DHX58(2), EIF4A3(2), ENOX1(2), G3BP1(2), HNRNPA2B1(2), HNRNPCL1(2), IGF2BP1(2), IGHMBP2(2), LIN28A(2), PABPC1(2), PABPC4(2), PPIE(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PIGQ in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
