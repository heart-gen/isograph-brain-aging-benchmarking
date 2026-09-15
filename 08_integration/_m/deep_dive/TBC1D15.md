# TBC1D15 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.502  ·  missense o/e 0.86

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | yes | yes | yes | no annotated structural change | risk allele T decreases usage of junction chr12:71918548-71920731(+) (ENST00000319106.12,E |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PPRC1, RBM14, DHX9, RBM41, SNRPB2, DDX58, ZC3H10, ADAR, RALY, RBM6, RBMS3, DHX58, CPEB2, IFIH1
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), ADAR(6), AGO1(6), AGO2(6), CELF4(6), CELF5(6), CNOT4(6), CPEB2(6), CPEB4(6), CSTF2(6), DDX19B(6), DDX58(6), DHX58(6), DHX9(6), EIF4A3(6), FMR1(6), G3BP2(6), GRSF1(6), HNRNPA1L2(6)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TBC1D15 in PD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
TBC1D15 is a Rab7 GTPase-activating protein at the mitochondria-lysosome interface, a pathway central to Parkinson's-disease mitophagy; its PD-associated splicing-led switch has no established isoform literature -- a mechanistically suggestive novel candidate.
