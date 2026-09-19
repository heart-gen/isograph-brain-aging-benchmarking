# IDH3B — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.10 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.037  ·  missense o/e 0.84

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | no | — | risk allele C decreases IDH3B expression (gene-level; no intron) |
| scz | eQTL | Brain_Frontal_Cortex_BA9 | 0.02 | no | — | yes | — | risk allele C decreases IDH3B expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | — | no | — | risk allele C increases usage of junction chr20:2658837-2659109(-) (unmapped transcript; n |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | no | — | risk allele C decreases usage of junction chr20:2658837-2659525(-) (ENST00000380843.9,ENST |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | no | — | risk allele C increases usage of junction chr20:2659312-2659525(-) (ENST00000613370.1; not |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | no | — | risk allele C decreases usage of junction chr20:2659793-2660030(-) (ENST00000380843.9,ENST |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | no | — | risk allele C increases usage of junction chr20:2659793-2660068(-) (ENST00000477689.2; not |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | — | yes | — | risk allele C increases usage of junction chr20:2658837-2659109(-) (unmapped transcript; n |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | yes | — | risk allele C decreases usage of junction chr20:2658837-2659525(-) (ENST00000380843.9,ENST |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | yes | — | risk allele C increases usage of junction chr20:2659312-2659525(-) (ENST00000613370.1; not |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | yes | — | risk allele C decreases usage of junction chr20:2659793-2660030(-) (ENST00000380843.9,ENST |
| scz | sQTL | Brain_Cerebellum | 0.10 | no | no | yes | — | risk allele C increases usage of junction chr20:2659793-2660068(-) (ENST00000477689.2; not |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Frontal_Cortex_BA9 | rs6115368 | C | -0.2718760073184967 | risk allele C decreases expression of IDH3B |
| scz | Brain_Cerebellum | rs6115368 | C | 1.5135231018066406 | risk allele C increases intron usage of chr20:2658837-2659109(-) |
| scz | Brain_Cerebellum | rs6115368 | C | -1.2718307971954346 | risk allele C decreases intron usage of chr20:2658837-2659525(-) |
| scz | Brain_Cerebellum | rs6115368 | C | 1.3039613962173462 | risk allele C increases intron usage of chr20:2659312-2659525(-) |
| scz | Brain_Cerebellum | rs6115368 | C | -0.641717791557312 | risk allele C decreases intron usage of chr20:2659793-2660030(-) |
| scz | Brain_Cerebellum | rs6115368 | C | 0.6258594989776611 | risk allele C increases intron usage of chr20:2659793-2660068(-) |
| scz | Brain_Frontal_Cortex_BA9 | rs6115368 | C | -0.2718760073184967 | risk allele C decreases expression of IDH3B |
| scz | Brain_Cerebellum | rs6115368 | C | 1.5135231018066406 | risk allele C increases intron usage of chr20:2658837-2659109(-) |
| scz | Brain_Cerebellum | rs6115368 | C | -1.2718307971954346 | risk allele C decreases intron usage of chr20:2658837-2659525(-) |
| scz | Brain_Cerebellum | rs6115368 | C | 1.3039613962173462 | risk allele C increases intron usage of chr20:2659312-2659525(-) |
| scz | Brain_Cerebellum | rs6115368 | C | -0.641717791557312 | risk allele C decreases intron usage of chr20:2659793-2660030(-) |
| scz | Brain_Cerebellum | rs6115368 | C | 0.6258594989776611 | risk allele C increases intron usage of chr20:2659793-2660068(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP1, G3BP2, ENOX1, MATR3, RBM46, YTHDC1, ZFP36L2, SAMD4A, SRSF5, DHX9, ESRP2, ADAR, RBM25, ESRP1, RBM5
- **All switched-motif RBPs (recurrence across regions):** AGO2(4), CPEB4(4), G3BP2(4), HNRNPC(4), G3BP1(4), ESRP1(4), PPIE(4), SRSF6(4), ZFP36L2(4), ZRANB2(4), HNRNPD(4), NONO(4), HNRNPU(4), HNRNPK(4), SRSF5(4), RBM46(4), RBM4(4), RBM25(4), YTHDC1(3), RBFOX1(3)

## 5. Interpretation
IDH3B has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
