# CTSH — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.39 (Brain_Hippocampus)
- **Constraint:** LOEUF 1.178  ·  missense o/e 0.91

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Substantia_nigra | 0.26 | no | — | yes | — | risk allele G increases CTSH expression (gene-level; no intron) |
| ad | sQTL | Brain_Hippocampus | 0.39 | no | no | yes | — | risk allele T decreases usage of junction chr15:78937423-78937686(-) (ENST00000533777.5,EN |
| ad | sQTL | Brain_Hippocampus | 0.39 | yes | yes | yes | no annotated structural change | risk allele T increases usage of junction chr15:78937423-78939140(-) (ENST00000220166.10,E |
| ad | sQTL | Brain_Hippocampus | 0.39 | no | no | yes | — | risk allele T decreases usage of junction chr15:78937823-78939140(-) (ENST00000533777.5,EN |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZCRB1, CPEB1, KHDRBS2, EIF4A3, DDX19B, U2AF2
- **All switched-motif RBPs (recurrence across regions):** CPEB1(4), CPEB4(4), DDX19B(4), EIF4A3(4), ZRANB2(4), ESRP2(4), HNRNPA3(4), IGF2BP1(4), PABPN1(4), KHDRBS2(4), RBFOX2(4), RBM28(4), ZCRB1(4), SNRPB2(4), U2AF2(4), ACO1(2), HNRNPLL(2), RBMY1A1(2), SAMD4A(2), YTHDC1(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for CTSH in AD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
