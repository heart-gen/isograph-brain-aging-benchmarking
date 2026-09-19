# THAP3 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.037  ·  missense o/e 0.90

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | yes | yes | no | no annotated structural change | risk allele T increases usage of junction chr1:6628691-6629643(+) (ENST00000487819.5; matc |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs3789572 | T | 0.5557003617286682 | risk allele T increases intron usage of chr1:6628691-6629643(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ACO1, SAMD4A, PABPN1, RBM25, CSTF2, SRSF10, HNRNPK, PABPC1, SART3, CNOT4, IGHMBP2, CPEB4, EIF4B, NOVA2, ZRANB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), CELF6(4), CPEB4(4), CPEB2(4), CPEB1(4), ELAVL4(4), EIF4B(4), CSTF2(4), HNRNPD(4), HNRNPC(4), HNRNPCL1(4), HNRNPAB(4), FXR1(4), KHDRBS2(4), KHDRBS3(4), KHDRBS1(4), IGHMBP2(4), IGF2BP2(4), HNRNPDL(4), IGF2BP3(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for THAP3 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.
