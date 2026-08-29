# PGS1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.09 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.998  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.09 | yes | yes | yes | no annotated structural change | risk allele A decreases usage of junction chr17:78378808-78392476(+) (ENST00000262764.11,E |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.09 | no | no | yes | — | risk allele A increases usage of junction chr17:78382785-78392476(+) (ENST00000586510.6; n |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM14, ZC3H10, ENOX1, CELF6, RALY, MATR3, ESRP1, ESRP2, HNRNPA3, A1CF, HNRNPCL1, AKAP1, SNRNP70, HNRNPA1L2
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), AGO2(4), AKAP1(4), CELF4(4), CELF5(4), CELF6(4), CPEB1(4), CPEB4(4), CSTF2(4), DDX19B(4), DHX58(4), EIF4A3(4), EIF4B(4), ELAVL3(4), ENOX1(4), ESRP1(4), ESRP2(4), G3BP1(4), GRSF1(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PGS1 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
PGS1 (phosphatidylglycerophosphate synthase 1; mitochondrial phospholipid biosynthesis) colocalizes as a splicing-led switch in ALS with no established disease isoform biology -- a novel candidate.
