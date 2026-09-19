# GSTO2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 1.238  ·  missense o/e 0.96

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | yes | yes | yes | no annotated structural change | risk allele A increases usage of junction chr10:104275334-104277894(+) (ENST00000338595.7, |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | yes | yes | yes | no annotated structural change | risk allele A decreases usage of junction chr10:104278116-104279370(+) (ENST00000338595.7, |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | yes | yes | yes | no annotated structural change | risk allele A increases usage of junction chr10:104278116-104297578(+) (ENST00000450629.6, |
| scz | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | yes | yes | yes | no annotated structural change | risk allele A decreases usage of junction chr10:104279471-104297578(+) (ENST00000338595.7, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs10883990 | A | 0.42206618189811707 | risk allele A increases intron usage of chr10:104275334-104277894(+) |
| scz | Brain_Frontal_Cortex_BA9 | rs10883990 | A | -0.841559112071991 | risk allele A decreases intron usage of chr10:104278116-104279370(+) |
| scz | Brain_Frontal_Cortex_BA9 | rs10883990 | A | 0.8467937707901001 | risk allele A increases intron usage of chr10:104278116-104297578(+) |
| scz | Brain_Frontal_Cortex_BA9 | rs10883990 | A | -0.8151730298995972 | risk allele A decreases intron usage of chr10:104279471-104297578(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPA0, KHDRBS1, HNRNPU, SYNCRIP, PABPC1, U2AF2, ELAVL3
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), ADAR(3), AGO1(3), AGO2(3), AKAP1(3), CELF6(3), CMTR1(3), CNOT4(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), DDX58(3), DHX58(3), DHX9(3), EIF4A3(3), EIF4B(3), ELAVL3(3), ENOX1(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for GSTO2 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
