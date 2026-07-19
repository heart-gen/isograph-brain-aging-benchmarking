# ARVCF — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.999  ·  missense o/e 1.14

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | no | yes | — | risk allele C decreases usage of junction chr22:19982091-19983673(-) (ENST00000462319.1; n |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr22:19982091-19990585(-) (ENST00000263207.8,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs4819527 | C | -0.5570746660232544 | risk allele C decreases intron usage of chr22:19982091-19983673(-) |
| scz | Brain_Cerebellar_Hemisphere | rs4819527 | C | 0.3532843291759491 | risk allele C increases intron usage of chr22:19982091-19990585(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PPRC1
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), CELF4(4), CELF5(4), CMTR1(4), DDX19B(4), DHX58(4), ELAVL3(4), ENOX1(4), ESRP2(4), FXR2(4), G3BP1(4), G3BP2(4), HNRNPA0(4), HNRNPAB(4), HNRNPD(4), HNRNPLL(4), IGHMBP2(4), KHDRBS1(4), KHDRBS2(4), KHDRBS3(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for ARVCF in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
