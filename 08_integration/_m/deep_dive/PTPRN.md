# PTPRN — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.735  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellum | 0.02 | yes | yes | yes | biotype switch, cds, coding status change, first | risk allele G decreases usage of junction chr2:219295141-219296226(-) (ENST00000295718.7,E |
| als | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele G increases usage of junction chr2:219295966-219296226(-) (unmapped transcript |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellum | rs6436132 | G | -0.6307845711708069 | risk allele G decreases intron usage of chr2:219295141-219296226(-) |
| als | Brain_Cerebellum | rs6436132 | G | 0.5719643235206604 | risk allele G increases intron usage of chr2:219295966-219296226(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, PPRC1, RBM8A, RC3H1, YBX2, ZC3H10, AKAP1, SRSF10, G3BP1, SNRPB2, CPEB1, TIA1, HNRNPA2B1
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), AGO1(5), AKAP1(5), AGO2(5), CELF6(5), CMTR1(5), CPEB4(5), CPEB1(5), CPEB2(5), ELAVL3(5), EIF4A3(5), EIF4B(5), DHX58(5), DDX19B(5), HNRNPC(5), HNRNPAB(5), HNRNPA0(5), HNRNPA3(5), HNRNPCL1(5), ENOX1(5)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for PTPRN in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
PTPRN/IA-2 is a dense-core vesicle transmembrane protein and a well-known autoantigen in type 1 diabetes, with a documented role in neuroendocrine secretion. Its isoform usage in brain, and in ALS, is not characterised.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
