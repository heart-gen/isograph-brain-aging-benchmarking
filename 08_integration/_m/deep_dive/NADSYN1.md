# NADSYN1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.929  ·  missense o/e 0.90

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele G increases usage of junction chr11:71448885-71451690(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele G decreases usage of junction chr11:71448885-71455110(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr11:71453381-71455110(+) (ENST00000319023.7,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr11:71474526-71478395(+) (ENST00000319023.7,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr11:71475622-71476201(+) (ENST00000528509.5; m |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr11:71477438-71478395(+) (ENST00000525200.5,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele G increases usage of junction chr11:71448885-71451690(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele G decreases usage of junction chr11:71448885-71455110(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr11:71453381-71455110(+) (ENST00000319023.7,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr11:71474526-71478395(+) (ENST00000319023.7,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr11:71475622-71476201(+) (ENST00000528509.5; m |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr11:71477438-71478395(+) (ENST00000525200.5,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | 0.8993797898292542 | risk allele G increases intron usage of chr11:71448885-71451690(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | -0.7316503524780273 | risk allele G decreases intron usage of chr11:71448885-71455110(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | -0.6603376269340515 | risk allele G decreases intron usage of chr11:71453381-71455110(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | -0.7384153604507446 | risk allele G decreases intron usage of chr11:71474526-71478395(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | 0.6088851094245911 | risk allele G increases intron usage of chr11:71475622-71476201(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | 0.6135239601135254 | risk allele G increases intron usage of chr11:71477438-71478395(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | 0.8993797898292542 | risk allele G increases intron usage of chr11:71448885-71451690(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | -0.7316503524780273 | risk allele G decreases intron usage of chr11:71448885-71455110(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | -0.6603376269340515 | risk allele G decreases intron usage of chr11:71453381-71455110(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | -0.7384153604507446 | risk allele G decreases intron usage of chr11:71474526-71478395(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | 0.6088851094245911 | risk allele G increases intron usage of chr11:71475622-71476201(+) |
| scz | Brain_Cerebellar_Hemisphere | rs12803256 | G | 0.6135239601135254 | risk allele G increases intron usage of chr11:71477438-71478395(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, RBM6, SNRPB2, RBMS3, HNRNPA3, CELF4, CELF5, FXR1, CELF6, MATR3, PABPC3, RBM46, PABPC5, YTHDC1, EIF4B
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), AGO1(4), AGO2(4), CMTR1(4), CELF6(4), CSTF2(4), DDX19B(4), CPEB4(4), CPEB2(4), ESRP1(4), EIF4B(4), EIF4A3(4), DDX58(4), YTHDC1(4), ZCRB1(4), ZRANB2(4), SUPV3L1(4), G3BP1(4), FXR1(4), HNRNPAB(4)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for NADSYN1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
NADSYN1 completes NAD+ synthesis and sits in a locus shared with DHCR7, which complicates gene attribution at the signal. Isoform biology is not established.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
