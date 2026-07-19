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
- **Switched *and* module-enriched (q<0.05) RBPs:** PABPC1, ESRP1, RNASEL, ELAVL4, HNRNPD, SYNCRIP, IGHMBP2, ELAVL3, SRSF10, DHX58, FXR2, SART3, KHDRBS3, U2AF2
- **All switched-motif RBPs (recurrence across regions):** DDX19B(3), DHX58(3), EIF4A3(3), ELAVL3(3), ELAVL4(3), ESRP1(3), FXR2(3), HNRNPA0(3), HNRNPD(3), IGHMBP2(3), KHDRBS2(3), KHDRBS3(3), PABPC1(3), PUM2(3), RNASEL(3), SART3(3), SNRPB2(3), SRSF10(3), SYNCRIP(3), U2AF2(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for RTEL1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
RTEL1 (telomere-maintenance helicase; AD/SCZ locus) has documented alternative C-terminal isoforms in other tissues, but a brain disease-specific splice role is not established -- a splicing-led candidate whose isoform choice warrants transcript-level follow-up.
