# IFNAR2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.761  ·  missense o/e 0.93

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellum | 0.03 | yes | yes | no | no annotated structural change | risk allele C increases usage of junction chr21:33245074-33246718(+) (ENST00000342101.7,EN |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | no | no | — | risk allele C decreases usage of junction chr21:33246482-33246718(+) (ENST00000682044.1,EN |
| ad | sQTL | Brain_Cerebellum | 0.03 | yes | yes | no | no annotated structural change | risk allele C decreases usage of junction chr21:33260727-33262793(+) (ENST00000342136.9,EN |
| ad | sQTL | Brain_Cerebellum | 0.03 | no | — | no | — | risk allele C increases usage of junction chr21:33276753-33279752(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellum | rs9975538 | C | 0.47209376096725464 | risk allele C increases intron usage of chr21:33245074-33246718(+) |
| ad | Brain_Cerebellum | rs9975538 | C | -0.38836008310317993 | risk allele C decreases intron usage of chr21:33246482-33246718(+) |
| ad | Brain_Cerebellum | rs9975538 | C | -0.6724862456321716 | risk allele C decreases intron usage of chr21:33260727-33262793(+) |
| ad | Brain_Cerebellum | rs9975538 | C | 0.5353890061378479 | risk allele C increases intron usage of chr21:33276753-33279752(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM42, ZFP36L2, RBM6, ZNF638, G3BP1, PPRC1, RBM8A, G3BP2, SNRPB2, CNOT4, RBM41, EIF4B, IGF2BP1, ZC3H10, AKAP1
- **All switched-motif RBPs (recurrence across regions):** ACO1(2), AGO2(2), AKAP1(2), CELF4(2), CELF5(2), CELF6(2), CNOT4(2), CPEB2(2), CPEB4(2), DDX19B(2), DHX58(2), EIF4A3(2), EIF4B(2), ENOX1(2), ESRP1(2), ESRP2(2), G3BP1(2), G3BP2(2), GRSF1(2), HNRNPCL1(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for IFNAR2 in AD. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.
