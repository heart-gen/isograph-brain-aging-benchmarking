# SPG7 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.433  ·  missense o/e 1.03

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | yes | no annotated structural change | risk allele G decreases usage of junction chr16:89544772-89545872(+) (ENST00000561945.2,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | yes | no annotated structural change | risk allele G increases usage of junction chr16:89544772-89546658(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | yes | no annotated structural change | risk allele G decreases usage of junction chr16:89548113-89550494(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G increases usage of junction chr16:89549713-89550494(+) (ENST00000569820.6,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G increases usage of junction chr16:89554563-89555878(+) (ENST00000645063.1,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | yes | no annotated structural change | risk allele G decreases usage of junction chr16:89554563-89556887(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | yes | no annotated structural change | risk allele G increases usage of junction chr16:89556085-89556887(+) (ENST00000565891.2,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | no | no annotated structural change | risk allele G decreases usage of junction chr16:89544772-89545872(+) (ENST00000561945.2,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | no | no annotated structural change | risk allele G increases usage of junction chr16:89544772-89546658(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | no | no annotated structural change | risk allele G decreases usage of junction chr16:89548113-89550494(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | no | — | risk allele G increases usage of junction chr16:89549713-89550494(+) (ENST00000569820.6,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | no | — | risk allele G increases usage of junction chr16:89554563-89555878(+) (ENST00000645063.1,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | no | no annotated structural change | risk allele G decreases usage of junction chr16:89554563-89556887(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | yes | yes | no | no annotated structural change | risk allele G increases usage of junction chr16:89556085-89556887(+) (ENST00000565891.2,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs12919215 | G | -0.45302918553352356 | risk allele G decreases intron usage of chr16:89544772-89545872(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 0.5540867447853088 | risk allele G increases intron usage of chr16:89544772-89546658(+) |
| scz | Brain_Cerebellum | rs12919215 | G | -0.9817131757736206 | risk allele G decreases intron usage of chr16:89548113-89550494(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 0.9417046308517456 | risk allele G increases intron usage of chr16:89549713-89550494(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 1.3150441646575928 | risk allele G increases intron usage of chr16:89554563-89555878(+) |
| scz | Brain_Cerebellum | rs12919215 | G | -1.6915165185928345 | risk allele G decreases intron usage of chr16:89554563-89556887(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 1.5796688795089722 | risk allele G increases intron usage of chr16:89556085-89556887(+) |
| scz | Brain_Cerebellum | rs12919215 | G | -0.45302918553352356 | risk allele G decreases intron usage of chr16:89544772-89545872(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 0.5540867447853088 | risk allele G increases intron usage of chr16:89544772-89546658(+) |
| scz | Brain_Cerebellum | rs12919215 | G | -0.9817131757736206 | risk allele G decreases intron usage of chr16:89548113-89550494(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 0.9417046308517456 | risk allele G increases intron usage of chr16:89549713-89550494(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 1.3150441646575928 | risk allele G increases intron usage of chr16:89554563-89555878(+) |
| scz | Brain_Cerebellum | rs12919215 | G | -1.6915165185928345 | risk allele G decreases intron usage of chr16:89554563-89556887(+) |
| scz | Brain_Cerebellum | rs12919215 | G | 1.5796688795089722 | risk allele G increases intron usage of chr16:89556085-89556887(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, G3BP1, RBM6, SNRPB2, G3BP2, ENOX1, ZC3H10, HNRNPA3, CELF4, CELF5, FXR1, CELF6, MATR3, SNRNP70
- **All switched-motif RBPs (recurrence across regions):** AGO2(6), DDX19B(6), CSTF2(6), CPEB2(6), G3BP2(6), HNRNPD(6), HNRNPA0(6), ESRP1(6), ELAVL3(6), EIF4B(6), EIF4A3(6), RBM41(6), RBM6(6), SNRPA(6), SNRNP70(6), SNRPB2(6), SRSF11(6), TUT1(6), SYNCRIP(6), PABPC4(6)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for SPG7 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
