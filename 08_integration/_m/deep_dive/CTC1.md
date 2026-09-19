# CTC1 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** ALS  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.05 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.776  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Caudate_basal_ganglia | 0.02 | no | — | no | — | risk allele G decreases CTC1 expression (gene-level; no intron) |
| als | sQTL | Brain_Cerebellum | 0.05 | no | no | no | — | risk allele G decreases usage of junction chr17:8228629-8229142(-) (ENST00000699853.1; not |
| als | sQTL | Brain_Cerebellum | 0.05 | yes | yes | no | no annotated structural change | risk allele G decreases usage of junction chr17:8230651-8231276(-) (ENST00000449476.7,ENST |
| als | sQTL | Brain_Cerebellum | 0.05 | yes | yes | no | no annotated structural change | risk allele G increases usage of junction chr17:8231469-8231726(-) (ENST00000449476.7,ENST |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Caudate_basal_ganglia | rs3027235 | G | -0.8031962513923645 | risk allele G decreases expression of CTC1 |
| als | Brain_Cerebellum | rs3027235 | G | -0.8849449753761292 | risk allele G decreases intron usage of chr17:8228629-8229142(-) |
| als | Brain_Cerebellum | rs3027235 | G | -2.087371587753296 | risk allele G decreases intron usage of chr17:8230651-8231276(-) |
| als | Brain_Cerebellum | rs3027235 | G | 2.0819549560546875 | risk allele G increases intron usage of chr17:8231469-8231726(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP2, ACO1, SAMD4A, ZFP36L2, ZNF638, SRSF11, CELF4, PABPC3, DHX9, CELF5, SNRPB2, PABPC1, SART3, ADAR, IGHMBP2
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), ADAR(2), AGO2(2), AKAP1(2), CELF4(2), CELF5(2), CELF6(2), CPEB2(2), CPEB4(2), DHX58(2), DHX9(2), ESRP1(2), G3BP1(2), YBX2(2), HNRNPA1L2(2), HNRNPA3(2), HNRNPAB(2), HNRNPLL(2), HNRNPDL(2)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for CTC1 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.
