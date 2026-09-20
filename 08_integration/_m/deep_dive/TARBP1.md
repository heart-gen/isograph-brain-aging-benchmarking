# TARBP1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.950  ·  missense o/e 1.12

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele T increases usage of junction chr1:234401262-234405312(-) (ENST00000468077.5;  |
| scz | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele T decreases usage of junction chr1:234401262-234405903(-) (ENST00000040877.2,E |
| scz | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele T increases usage of junction chr1:234405480-234405903(-) (ENST00000468077.5;  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs1039996 | T | 1.312606692314148 | risk allele T increases intron usage of chr1:234401262-234405312(-) |
| scz | Brain_Cerebellum | rs1039996 | T | -1.449627161026001 | risk allele T decreases intron usage of chr1:234401262-234405903(-) |
| scz | Brain_Cerebellum | rs1039996 | T | 1.312963604927063 | risk allele T increases intron usage of chr1:234405480-234405903(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** DDX58, DHX9, PABPC5, ZNF346, PPRC1, ADAR, IFIH1, ERI1, CNOT4, RBM46, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), ADAR(2), AGO1(2), AGO2(2), AKAP1(2), CELF4(2), CELF5(2), CMTR1(2), CNOT4(2), CPEB2(2), CPEB4(2), CSTF2(2), DDX19B(2), DDX58(2), DHX58(2), DHX9(2), EIF4A3(2), EIF4B(2), ERI1(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for TARBP1 in SCZ. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
TARBP1 is a TRBP-related RNA methyltransferase with no established disease isoform biology; its usage range touches the detection floor here.

_Curation: novel_candidate._ No established disease-specific isoform biology was found. That is a statement about the literature, not about the evidence here: an under-characterised switch is what this method is built to surface, and this row is a nomination rather than a null result.
