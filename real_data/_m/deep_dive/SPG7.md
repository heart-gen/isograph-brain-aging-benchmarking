# SPG7 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.04 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.433  ·  missense o/e 1.03

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G decreases usage of junction chr16:89544772-89545872(+) (ENST00000561945.2,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G increases usage of junction chr16:89544772-89546658(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G decreases usage of junction chr16:89548113-89550494(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G increases usage of junction chr16:89549713-89550494(+) (ENST00000569820.6,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G increases usage of junction chr16:89554563-89555878(+) (ENST00000645063.1,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G decreases usage of junction chr16:89554563-89556887(+) (ENST00000268704.7,EN |
| scz | sQTL | Brain_Cerebellum | 0.04 | no | no | yes | — | risk allele G increases usage of junction chr16:89556085-89556887(+) (ENST00000565891.2,EN |

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

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SNRPB2, HNRNPCL1, SART3, CPEB2, RALY, G3BP2, YBX2, PTBP2, ENOX1, ZFP36L2
- **All switched-motif RBPs (recurrence across regions):** CELF5(3), CELF4(3), CSTF2(3), CPEB2(3), DDX19B(3), ELAVL3(3), EIF4B(3), EIF4A3(3), ZFP36L2(3), ZC3H10(3), RBM6(3), ESRP1(3), HNRNPAB(3), HNRNPCL1(3), HNRNPA0(3), G3BP2(3), KHDRBS2(3), SNRPB2(3), SRSF11(3), YBX2(3)

## 5. Interpretation
SPG7 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
