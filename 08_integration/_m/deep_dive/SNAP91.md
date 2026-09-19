# SNAP91 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.05 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.362  ·  missense o/e 0.96

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.05 | no | no | yes | — | risk allele T increases usage of junction chr6:83593017-83593478(-) (ENST00000518312.5; no |
| scz | sQTL | Brain_Cerebellum | 0.05 | yes | yes | yes | no annotated structural change | risk allele T increases usage of junction chr6:83607808-83610650(-) (ENST00000195649.10,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs217342 | T | 0.42374029755592346 | risk allele T increases intron usage of chr6:83593017-83593478(-) |
| scz | Brain_Cerebellum | rs217342 | T | 0.7560865879058838 | risk allele T increases intron usage of chr6:83607808-83610650(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, PPRC1, RBM8A, SYNCRIP, G3BP1
- **All switched-motif RBPs (recurrence across regions):** DHX58(7), FXR2(7), G3BP1(7), HNRNPLL(7), IGHMBP2(7), LIN28A(7), MATR3(7), PABPC4(7), PABPN1(7), PPRC1(7), RBM14(7), RBM28(7), RBM8A(7), SYNCRIP(7)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for SNAP91 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
