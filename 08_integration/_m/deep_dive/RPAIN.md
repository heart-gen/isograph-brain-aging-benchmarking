# RPAIN — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cortex)
- **Constraint:** LOEUF 1.320  ·  missense o/e 0.89

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cortex | 0.02 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr17:5426299-5428071(+) (ENST00000381209.8,ENST |
| scz | sQTL | Brain_Cortex | 0.02 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr17:5428211-5432542(+) (ENST00000381209.8,ENST |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs2189336 | G | -0.3790275454521179 | risk allele G decreases intron usage of chr17:5426299-5428071(+) |
| scz | Brain_Cortex | rs2189336 | G | -0.36680907011032104 | risk allele G decreases intron usage of chr17:5428211-5432542(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** NOVA1 , KHDRBS1, AGO1, RBMY1A1, NOVA2, PABPC1, TRA2A
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), AGO1(2), AKAP1(2), CPEB2(2), CELF6(2), HNRNPM(2), HNRNPA3(2), CPEB4(2), CSTF2(2), DHX58(2), DDX19B(2), EIF4A3(2), FXR2(2), HNRNPLL(2), HNRNPCL1(2), HNRNPAB(2), NELFE(2), KHDRBS3(2), MATR3(2), KHDRBS2(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for RPAIN in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.

## 6. Literature (known isoform biology)
RPAIN/RIP is an RPA-interacting nuclear import factor with no established brain disease isoform biology. Its junction validates in the PSI catalogue and in the recount, so the switch is well measured and the literature is simply absent.

_Curation: novel_candidate._ No established disease-specific isoform biology was found. That is a statement about the literature, not about the evidence here: an under-characterised switch is what this method is built to surface, and this row is a nomination rather than a null result.
