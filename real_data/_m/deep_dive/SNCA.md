# SNCA — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** LBD,PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cortex)
- **Constraint:** LOEUF 0.397  ·  missense o/e 0.71

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | sQTL | Brain_Cortex | 0.04 | yes | yes | yes | no annotated structural change | risk allele A increases usage of junction chr4:89835692-89836127(-) (ENST00000508895.5,ENS |
| pd | sQTL | Brain_Frontal_Cortex_BA9 | 0.03 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr4:89835692-89836127(-) (ENST00000508895.5,ENS |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB2, A1CF, DHX9, RBMS3, ZCRB1, ZFP36L2, RBMS1, AKAP1, RALY, ADAR, SUPV3L1, PABPC4, PABPC5, RBM41, RBM42
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ADAR(5), AGO1(5), AGO2(5), AKAP1(5), CELF4(5), CPEB2(5), EIF4B(5), DHX9(5), RBMS3(5), SUPV3L1(5), ZRANB2(5), PABPC4(5), PHAX(5), PTBP2(5), PUM2(5), PABPC5(5), LIN28A(5), IFIH1(5), ZCRB1(5)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for SNCA in LBD,PD — the same switch is genetically anchored across more than one trait. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
