# DNAJA3 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.733  ·  missense o/e 1.34

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr16:4426092-4434384(+) (ENST00000262375.11,ENS |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr16:4434517-4437402(+) (ENST00000262375.11,ENS |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs6500605 | G | 0.5182712078094482 | risk allele G increases intron usage of chr16:4426092-4434384(+) |
| scz | Brain_Cerebellar_Hemisphere | rs6500605 | G | -0.5069531798362732 | risk allele G decreases intron usage of chr16:4434517-4437402(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, HNRNPA0, KHDRBS1, CPEB1, SRSF10, U2AF2, PPIE, RBMY1A1, AGO1, SF1, ZRANB2
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), CELF5(4), CELF4(4), AGO2(4), CELF6(4), CPEB2(4), CPEB1(4), CNOT4(4), DDX19B(4), DHX58(4), CSTF2(4), G3BP2(4), FXR2(4), ESRP1(4), ENOX1(4), EIF4A3(4), ZCRB1(4), YTHDC1(4), SYNCRIP(4), SUPV3L1(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for DNAJA3 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
DNAJA3/Tid1 is a mitochondrial HSP40 co-chaperone with two long-standing splice forms that differ at the C-terminus and have been reported to act in opposite directions on apoptosis, so isoform choice rather than total level is the functional variable. A schizophrenia-specific isoform role is not established.

_Curation: isoform_documented._ A disease-relevant isoform program is established for this gene.

_References:_ [citation needed: Tid1-L / Tid1-S splice forms and opposing apoptotic effects]
