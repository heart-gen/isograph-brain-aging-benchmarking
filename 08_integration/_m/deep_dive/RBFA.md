# RBFA — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.977  ·  missense o/e 1.03

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr18:80042219-80044212(+) (ENST00000306735.10,E |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr18:80048193-80069798(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele C decreases usage of junction chr18:80067308-80069798(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs8098075 | C | 0.5398757457733154 | risk allele C increases intron usage of chr18:80042219-80044212(+) |
| scz | Brain_Cerebellum | rs8098075 | C | 0.6476700305938721 | risk allele C increases intron usage of chr18:80048193-80069798(+) |
| scz | Brain_Cerebellum | rs8098075 | C | -1.2658482789993286 | risk allele C decreases intron usage of chr18:80067308-80069798(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, RBM41, SNRPB2, RBMS3, PPRC1, ENOX1, IGF2BP1, HNRNPA3, CELF4, CELF5, CELF6, PABPC3, PABPC5, YTHDC1, EIF4B
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO1(3), AGO2(3), AKAP1(3), ANKHD1(3), CELF4(3), CELF5(3), CELF6(3), CMTR1(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), DHX58(3), EIF4A3(3), EIF4B(3), ELAVL3(3), ENOX1(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for RBFA in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.
