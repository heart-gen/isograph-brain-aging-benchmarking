# TBC1D15 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.502  ·  missense o/e 0.86

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele T decreases usage of junction chr12:71918548-71920731(+) (ENST00000319106.12,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Frontal_Cortex_BA9 | rs61754230 | T | -1.7267097234725952 | risk allele T decreases intron usage of chr12:71918548-71920731(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, RBM6, SNRPB2, RBMS3, PPRC1, G3BP2, CNOT4, ZC3H10, HNRNPA3, CELF4, CELF5, MATR3, PABPC3, RBM46
- **All switched-motif RBPs (recurrence across regions):** A1CF(7), ACO1(7), ADAR(7), AGO1(7), AGO2(7), ANKHD1(7), CELF4(7), CELF5(7), CNOT4(7), CPEB2(7), CPEB4(7), CSTF2(7), DDX19B(7), DDX58(7), DHX58(7), DHX9(7), EIF4A3(7), FMR1(7), G3BP2(7), GRSF1(7)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TBC1D15 in PD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.

## 6. Literature (known isoform biology)
TBC1D15 is a Rab7 GTPase-activating protein at the mitochondria-lysosome interface, a pathway central to Parkinson's mitophagy. Its PD-associated splicing-led switch has no established isoform literature -- a mechanistically suggestive nomination rather than a confirmation.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
