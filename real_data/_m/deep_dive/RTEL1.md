# RTEL1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** AD,SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.638  ·  missense o/e 1.10

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele C increases usage of junction chr20:63661496-63661762(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | yes | no annotated structural change | risk allele C decreases usage of junction chr20:63690947-63691742(+) (ENST00000318100.9,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs190464636 | C | -1.9781608581542969 | risk allele C decreases intron usage of chr20:63690947-63691742(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** IGF2BP1, RBM6, ESRP1, RNASEL, YTHDC1, PABPC1, RBM24, ELAVL4, DHX58, HNRNPC, HNRNPCL1, RALY
- **All switched-motif RBPs (recurrence across regions):** CPEB4(4), EIF4A3(4), DHX58(4), DDX19B(4), CSTF2(4), ELAVL3(4), ELAVL4(4), ESRP1(4), FXR2(4), HNRNPD(4), HNRNPCL1(4), HNRNPC(4), HNRNPA0(4), YTHDC1(4), ZFP36L2(4), U2AF2(4), TRA2A(4), HNRNPU(4), KHDRBS3(4), KHDRBS2(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for RTEL1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
RTEL1 (telomere-maintenance helicase; AD/SCZ locus) has documented alternative C-terminal isoforms in other tissues, but a brain disease-specific splice role is not established -- a splicing-led candidate whose isoform choice warrants transcript-level follow-up.
