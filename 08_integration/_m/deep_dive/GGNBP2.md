# GGNBP2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.06 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.201  ·  missense o/e 0.66

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.06 | no | — | no | — | risk allele C increases GGNBP2 expression (gene-level; no intron) |
| als | sQTL | Brain_Cerebellum | 0.01 | no | no | no | — | risk allele C increases usage of junction chr17:36560871-36567663(+) (ENST00000613102.5,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Spinal_cord_cervical_c-1 | rs9903355 | C | 0.4017837345600128 | risk allele C increases expression of GGNBP2 |
| als | Brain_Cerebellum | rs6607326 | C | 0.42049580812454224 | risk allele C increases intron usage of chr17:36560871-36567663(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, G3BP1, SNRPB2, RBMS3, G3BP2, ENOX1, CELF4, CELF5, FXR1, CELF6, MATR3, PABPC3, RBM46, PABPC5, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), AGO2(4), AKAP1(4), CELF5(4), CELF4(4), CELF6(4), CPEB2(4), DHX58(4), CSTF2(4), DDX19B(4), DDX58(4), ENOX1(4), EIF4A3(4), FXR1(4), ESRP1(4), U2AF2(4), YTHDC1(4), FXR2(4), G3BP1(4), GRSF1(4)

## 5. Interpretation
GGNBP2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.

## 6. Literature (known isoform biology)
GGNBP2/ZNF403 (17q12) is LoF-constrained (LOEUF 0.20) and colocalizes as a splicing-led switch in ALS; its isoform biology in neurodegeneration is uncharacterized -- a novel splicing-led candidate.
