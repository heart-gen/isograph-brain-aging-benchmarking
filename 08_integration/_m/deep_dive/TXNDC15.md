# TXNDC15 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.896  ·  missense o/e 0.79

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Caudate_basal_ganglia | 0.03 | no | no | no | — | risk allele T decreases usage of junction chr5:134874530-134893492(+) (ENST00000511070.5;  |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Caudate_basal_ganglia | rs58746189 | T | -0.7040115594863892 | risk allele T decreases intron usage of chr5:134874530-134893492(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SAMD4A, A1CF, IGF2BP1, RBM6, G3BP1, RBM41, SNRPB2, RBM25, CELF4, DHX9, CSTF2, SRSF10, YTHDC1, CELF5, ADAR
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ADAR(3), AGO2(3), CELF4(3), CELF5(3), CELF6(3), CMTR1(3), CPEB2(3), CPEB4(3), CSTF2(3), DDX19B(3), DDX58(3), DHX9(3), EIF4A3(3), EIF4B(3), ELAVL3(3), ENOX1(3), FXR2(3), G3BP1(3), GRSF1(3)

## 5. Interpretation
TXNDC15 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
