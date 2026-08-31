# ARVCF — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.999  ·  missense o/e 1.14

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | yes | — | risk allele C decreases usage of junction chr22:19982091-19983673(-) (ENST00000462319.1; n |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr22:19982091-19990585(-) (ENST00000263207.8,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs4819527 | C | -0.5570746660232544 | risk allele C decreases intron usage of chr22:19982091-19983673(-) |
| scz | Brain_Cerebellar_Hemisphere | rs4819527 | C | 0.3532843291759491 | risk allele C increases intron usage of chr22:19982091-19990585(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** KHDRBS1, HNRNPD, ELAVL3, PABPC1, RNASEL, PPIE, ZFP36, CPEB4, TIAL1, HNRNPDL, AGO1, NUDT21, U2AF2
- **All switched-motif RBPs (recurrence across regions):** ACO1(6), AGO1(6), AGO2(6), CELF4(6), CELF5(6), CMTR1(6), CNOT4(6), CPEB2(6), CPEB4(6), CSTF2(6), DAZAP1(6), DDX19B(6), DHX58(6), EIF4A3(6), EIF4B(6), ELAVL2(6), ELAVL3(6), ENOX1(6), ESRP2(6), FXR1(6)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for ARVCF in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
ARVCF lies in the 22q11.2 schizophrenia deletion region (haplotypic SCZ association with COMT) and is itself a modulator of pre-mRNA splicing -- it interacts with SRSF1, DDX5 and hnRNP H2 and alters alternative-splicing activity -- so a splicing-led ARVCF event is consistent with both its locus and its molecular role.

_References:_ @doi:10.1038/sj.mp.4001586
