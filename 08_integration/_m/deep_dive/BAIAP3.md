# BAIAP3 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Hypothalamus)
- **Constraint:** LOEUF 1.383  ·  missense o/e 1.17

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Hypothalamus | 0.01 | yes | yes | no | biotype switch, cds, coding status change, first | risk allele G increases usage of junction chr16:1333749-1338540(+) (ENST00000397488.6,ENST |
| als | sQTL | Brain_Hypothalamus | 0.01 | no | no | no | — | risk allele G decreases usage of junction chr16:1334756-1338540(+) (ENST00000324385.9,ENST |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Hypothalamus | rs185060782 | G | 2.091409921646118 | risk allele G increases intron usage of chr16:1333749-1338540(+) |
| als | Brain_Hypothalamus | rs185060782 | G | -1.5065308809280396 | risk allele G decreases intron usage of chr16:1334756-1338540(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PPRC1, HNRNPLL, RBM14, PABPC3, YTHDC1, RBM8A
- **All switched-motif RBPs (recurrence across regions):** HNRNPA3(4), HNRNPLL(4), HNRNPM(4), PABPC3(4), RBM14(4), PPRC1(4), RBM4(4), RBM8A(4), ZC3H10(4), SNRPA(4), SRSF10(4), SRSF4(4), ZNF638(4), YTHDC1(4), ENOX1(1), SNRNP70(1)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for BAIAP3 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage.

## 6. Literature (known isoform biology)
BAIAP3 is a Munc13-family protein controlling dense-core vesicle secretion, so it sits in the same secretory machinery as several other genes in this panel. Isoform biology is not established, and its usage range here touches the detection floor.

_Curation: gene_documented._ The gene and its disease association are established; which isoform the risk variant selects is not characterised.
