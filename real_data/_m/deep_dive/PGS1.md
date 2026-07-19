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
- **Switched *and* module-enriched (q<0.05) RBPs:** RALY, HNRNPCL1, CELF6, ZFP36L2, ENOX1, HNRNPA1L2, SNRNP70, IGF2BP2, ESRP1, AKAP1, ESRP2, DHX58, EIF4A3, CELF5, CELF4
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO2(3), AKAP1(3), CELF4(3), CELF5(3), CELF6(3), CPEB1(3), CPEB4(3), CSTF2(3), DDX19B(3), DHX58(3), EIF4A3(3), ELAVL3(3), ENOX1(3), ESRP1(3), ESRP2(3), G3BP1(3), HNRNPA0(3), HNRNPA1L2(3), HNRNPA3(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PGS1 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
PGS1 (phosphatidylglycerophosphate synthase 1; mitochondrial phospholipid biosynthesis) colocalizes as a splicing-led switch in ALS with no established disease isoform biology -- a novel candidate.
