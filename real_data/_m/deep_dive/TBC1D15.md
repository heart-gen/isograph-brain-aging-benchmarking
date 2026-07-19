# TBC1D15 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Frontal_Cortex_BA9)
- **Constraint:** LOEUF 0.502  ·  missense o/e 0.86

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Frontal_Cortex_BA9 | 0.02 | yes | yes | yes | no annotated structural change | risk allele T decreases usage of junction chr12:71918548-71920731(+) (ENST00000319106.12,E |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB2, A1CF, DHX9, RBMS3, ZFP36L2, RBM41, RALY, ADAR, DHX58, PABPC5, PABPC4, MATR3, DDX58, SART3, RBM46
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), ADAR(4), AGO1(4), AGO2(4), CELF4(4), CELF5(4), CNOT4(4), CPEB2(4), CPEB4(4), CSTF2(4), DDX58(4), DHX58(4), DHX9(4), EIF4A3(4), EIF4B(4), G3BP2(4), HNRNPA1L2(4), HNRNPAB(4), HNRNPCL1(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TBC1D15 in PD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
TBC1D15 is a Rab7 GTPase-activating protein at the mitochondria-lysosome interface, a pathway central to Parkinson's-disease mitophagy; its PD-associated splicing-led switch has no established isoform literature -- a mechanistically suggestive novel candidate.
