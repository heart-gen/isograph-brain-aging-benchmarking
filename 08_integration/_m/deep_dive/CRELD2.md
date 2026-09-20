# CRELD2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.06 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.090  ·  missense o/e 1.00

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.06 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr22:49921761-49922294(+) (ENST00000404488.7; m |
| scz | sQTL | Brain_Cerebellum | 0.06 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr22:49921761-49922297(+) (ENST00000404488.7; m |
| scz | sQTL | Brain_Cerebellum | 0.06 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr22:49921761-49922612(+) (ENST00000328268.9,EN |
| scz | sQTL | Brain_Cerebellum | 0.06 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr22:49922443-49922612(+) (ENST00000404488.7; m |
| scz | sQTL | Brain_Cerebellum | 0.06 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr22:49924455-49925417(+) (ENST00000328268.9,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs7284417 | G | -1.037359356880188 | risk allele G decreases intron usage of chr22:49921761-49922294(+) |
| scz | Brain_Cerebellum | rs7284417 | G | -0.8943570256233215 | risk allele G decreases intron usage of chr22:49921761-49922297(+) |
| scz | Brain_Cerebellum | rs7284417 | G | 1.2438035011291504 | risk allele G increases intron usage of chr22:49921761-49922612(+) |
| scz | Brain_Cerebellum | rs7284417 | G | -1.1440813541412354 | risk allele G decreases intron usage of chr22:49922443-49922612(+) |
| scz | Brain_Cerebellum | rs7284417 | G | 0.5303111672401428 | risk allele G increases intron usage of chr22:49924455-49925417(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPA0, ELAVL3, KHDRBS1, HNRNPD, IGF2BP3, CPEB1, SRSF10, PABPC1, CPEB4, U2AF2, PPIE, RBMY1A1, PUM1, NOVA2, OAS1
- **All switched-motif RBPs (recurrence across regions):** AGO2(3), CELF6(3), CNOT4(3), CPEB1(3), CPEB4(3), CSTF2(3), DHX58(3), ELAVL3(3), ENOX1(3), ESRP1(3), ESRP2(3), G3BP2(3), HNRNPA0(3), HNRNPA3(3), HNRNPAB(3), HNRNPD(3), HNRNPK(3), HNRNPU(3), IGF2BP3(3), IGHMBP2(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for CRELD2 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
CRELD2 is an ER-stress-responsive secreted protein induced through the ATF6 arm of the unfolded protein response. It has the broadest SMR support in the panel (48 of 51 probes) and no documented isoform program.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
