# ADAM15 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** LBD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Cortex)
- **Constraint:** LOEUF 0.750  ·  missense o/e 0.86

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| lbd | sQTL | Brain_Cortex | 0.03 | no | no | no | — | risk allele A decreases usage of junction chr1:155060832-155061415(+) (ENST00000355956.6,E |
| lbd | sQTL | Brain_Cortex | 0.03 | no | no | no | — | risk allele A decreases usage of junction chr1:155060832-155061418(+) (ENST00000449910.6,E |
| lbd | sQTL | Brain_Cortex | 0.03 | no | no | no | — | risk allele A decreases usage of junction chr1:155060832-155062245(+) (ENST00000271836.10, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| lbd | Brain_Cortex | rs11589479 | A | -0.7571733593940735 | risk allele A decreases intron usage of chr1:155060832-155061415(+) |
| lbd | Brain_Cortex | rs11589479 | A | -1.100714087486267 | risk allele A decreases intron usage of chr1:155060832-155061418(+) |
| lbd | Brain_Cortex | rs11589479 | A | -0.650272786617279 | risk allele A decreases intron usage of chr1:155060832-155062245(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM6, RBMS3, ENOX1, IGF2BP1, CELF4, CELF5, CELF6, MATR3, SNRNP70, YTHDC1, ACO1, SAMD4A, RBM24, ESRP2
- **All switched-motif RBPs (recurrence across regions):** CPEB4(6), CELF6(6), MATR3(6), HNRNPU(6), HNRNPDL(6), NELFE(5), HNRNPAB(5), IGF2BP1(5), ESRP2(5), HNRNPC(5), DAZAP1(5), RBM6(5), ZRANB2(5), IGF2BP3(4), ELAVL4(4), ELAVL2(4), DHX58(4), AGO2(4), ENOX1(4), CELF5(4)

## 5. Interpretation
ADAM15 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
