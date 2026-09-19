# FLCN — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** AD,SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.487  ·  missense o/e 0.81

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr17:17213856-17214985(-) (ENST00000285071.9; m |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C decreases usage of junction chr17:17215316-17217069(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr17:17221459-17221537(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr17:17224143-17224458(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | yes | no annotated structural change | risk allele C decreases usage of junction chr17:17224143-17226176(-) (ENST00000285071.9,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr17:17225207-17226176(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | yes | yes | yes | no annotated structural change | risk allele C decreases usage of junction chr17:17228161-17231794(-) (ENST00000285071.9,EN |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr17:17228821-17229488(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr17:17228821-17231794(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C decreases usage of junction chr17:17228905-17231794(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr17:17229708-17231794(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T increases usage of junction chr17:17207069-17207539(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | no annotated structural change | risk allele T decreases usage of junction chr17:17213856-17214985(-) (ENST00000285071.9; m |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T increases usage of junction chr17:17215316-17217069(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17221459-17221537(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17224143-17224458(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | no annotated structural change | risk allele T increases usage of junction chr17:17224143-17226176(-) (ENST00000285071.9,EN |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17224511-17226176(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17225207-17226176(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | no annotated structural change | risk allele T increases usage of junction chr17:17228161-17231794(-) (ENST00000285071.9,EN |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17228821-17229488(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17228821-17231794(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T increases usage of junction chr17:17228905-17231794(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17229708-17231794(-) (unmapped transcript; |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T decreases usage of junction chr17:17231635-17231794(-) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellum | rs1708618 | T | 0.5260412693023682 | risk allele T increases intron usage of chr17:17207069-17207539(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.34788498282432556 | risk allele T decreases intron usage of chr17:17213856-17214985(-) |
| ad | Brain_Cerebellum | rs1708618 | T | 0.2575632631778717 | risk allele T increases intron usage of chr17:17215316-17217069(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.4006943106651306 | risk allele T decreases intron usage of chr17:17221459-17221537(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.543351411819458 | risk allele T decreases intron usage of chr17:17224143-17224458(-) |
| ad | Brain_Cerebellum | rs1708618 | T | 0.5833153128623962 | risk allele T increases intron usage of chr17:17224143-17226176(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.48275724053382874 | risk allele T decreases intron usage of chr17:17224511-17226176(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.6422919034957886 | risk allele T decreases intron usage of chr17:17225207-17226176(-) |
| ad | Brain_Cerebellum | rs1708618 | T | 0.7682068943977356 | risk allele T increases intron usage of chr17:17228161-17231794(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.8981060981750488 | risk allele T decreases intron usage of chr17:17228821-17229488(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.905381441116333 | risk allele T decreases intron usage of chr17:17228821-17231794(-) |
| ad | Brain_Cerebellum | rs1708618 | T | 0.6389769315719604 | risk allele T increases intron usage of chr17:17228905-17231794(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -1.0702159404754639 | risk allele T decreases intron usage of chr17:17229708-17231794(-) |
| ad | Brain_Cerebellum | rs1708618 | T | -0.36361023783683777 | risk allele T decreases intron usage of chr17:17231635-17231794(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | 0.37689974904060364 | risk allele C increases intron usage of chr17:17213856-17214985(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | -0.167857825756073 | risk allele C decreases intron usage of chr17:17215316-17217069(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | 0.30683213472366333 | risk allele C increases intron usage of chr17:17221459-17221537(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | 0.6497754454612732 | risk allele C increases intron usage of chr17:17224143-17224458(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | -0.5650355219841003 | risk allele C decreases intron usage of chr17:17224143-17226176(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | 0.40406084060668945 | risk allele C increases intron usage of chr17:17225207-17226176(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | -0.5289745926856995 | risk allele C decreases intron usage of chr17:17228161-17231794(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | 0.7103748917579651 | risk allele C increases intron usage of chr17:17228821-17229488(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | 0.8939891457557678 | risk allele C increases intron usage of chr17:17228821-17231794(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | -0.6771350502967834 | risk allele C decreases intron usage of chr17:17228905-17231794(-) |
| scz | Brain_Cerebellar_Hemisphere | rs1708618 | C | 0.9634525179862976 | risk allele C increases intron usage of chr17:17229708-17231794(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, G3BP1, SNRPB2, RBMS3, G3BP2, ENOX1, ZC3H10, HNRNPA3, CELF4, CELF5, FXR1, CELF6, MATR3, SNRNP70, RBM46
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), AGO1(5), AGO2(5), CPEB1(5), CELF4(5), CELF5(5), CELF6(5), CPEB4(5), CPEB2(5), CSTF2(5), DAZAP1(5), G3BP2(5), DDX19B(5), EIF4B(5), EIF4A3(5), ELAVL3(5), ENOX1(5), ESRP2(5), ESRP1(5)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for FLCN in AD,SCZ — the same switch is genetically anchored across more than one trait. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
