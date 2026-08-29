# NDRG4 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Anterior_cingulate_cortex_BA24)
- **Constraint:** LOEUF 0.756  ·  missense o/e 0.93

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.03 | no | no | yes | — | risk allele T decreases usage of junction chr16:58509352-58510645(+) (ENST00000421602.6,EN |
| scz | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.03 | no | no | yes | — | risk allele T increases usage of junction chr16:58509352-58511422(+) (ENST00000258187.9,EN |
| scz | sQTL | Brain_Anterior_cingulate_cortex_BA24 | 0.03 | no | no | yes | — | risk allele T decreases usage of junction chr16:58510683-58511422(+) (ENST00000421602.6,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs149165 | T | -0.33587732911109924 | risk allele T decreases intron usage of chr16:58509352-58510645(+) |
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs149165 | T | 0.34704241156578064 | risk allele T increases intron usage of chr16:58509352-58511422(+) |
| scz | Brain_Anterior_cingulate_cortex_BA24 | rs149165 | T | -0.3566074073314667 | risk allele T decreases intron usage of chr16:58510683-58511422(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CPEB1
- **All switched-motif RBPs (recurrence across regions):** CELF6(4), DHX58(4), IGF2BP2(4), HNRNPA1L2(4), FXR2(4), ESRP2(4), ENOX1(4), LIN28A(4), SRSF11(4), YBX2(4), YTHDC1(4), RBMY1A1(4), IGF2BP1(4), HNRNPM(4), HNRNPCL1(4), RBM25(4), RBM4(4), PABPC5(4), CPEB1(2), ACO1(2)

## 5. Interpretation
NDRG4 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
