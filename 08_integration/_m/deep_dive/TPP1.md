# TPP1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.48 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.746  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellum | 0.48 | no | — | yes | — | risk allele G decreases usage of junction chr11:6611818-6611958(-) (unmapped transcript; n |
| als | sQTL | Brain_Cerebellum | 0.48 | yes | yes | yes | no annotated structural change | risk allele G increases usage of junction chr11:6614686-6614866(-) (ENST00000299427.12,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellum | rs2072651 | G | -0.5992392301559448 | risk allele G decreases intron usage of chr11:6611818-6611958(-) |
| als | Brain_Cerebellum | rs2072651 | G | 0.6017784476280212 | risk allele G increases intron usage of chr11:6614686-6614866(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CSTF2, SRSF10, YTHDC1, SNRPB2, SART3, MATR3, RBM6, ZC3H10, PUM2, SRSF4, AKAP1, RBMS3, G3BP1, CELF6
- **All switched-motif RBPs (recurrence across regions):** AKAP1(4), CELF6(4), HNRNPLL(4), ESRP1(4), HNRNPAB(4), HNRNPA1L2(4), SRSF4(4), TARDBP(4), ZC3H10(4), RBMS3(4), NELFE(4), HNRNPM(4), ZCRB1(4), YTHDC1(4), SNRPB2(4), SART3(4), PABPC5(3), AGO2(3), RBM6(3), CPEB1(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TPP1 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
