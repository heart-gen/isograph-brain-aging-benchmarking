# PLXNB2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.545  ·  missense o/e 0.81

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | no | — | risk allele C decreases usage of junction chr22:50275963-50276629(-) (ENST00000359337.9,EN |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | no | — | risk allele C increases usage of junction chr22:50281499-50281566(-) (ENST00000359337.9,EN |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | — | no | — | risk allele C decreases usage of junction chr22:50281499-50281643(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | — | no | — | risk allele C decreases usage of junction chr22:50290597-50307553(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | — | no | — | risk allele C increases usage of junction chr22:50294778-50295377(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | no | — | risk allele C decreases usage of junction chr22:50294778-50301347(-) (ENST00000432455.5; n |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | no | — | risk allele C decreases usage of junction chr22:50294778-50307553(-) (ENST00000359337.9,EN |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | — | no | — | risk allele C increases usage of junction chr22:50295429-50301347(-) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | — | no | — | risk allele C increases usage of junction chr22:50295429-50307553(-) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs111925803 | C | -0.40249067544937134 | risk allele C decreases intron usage of chr22:50275963-50276629(-) |
| scz | Brain_Cerebellum | rs111925803 | C | 0.48778530955314636 | risk allele C increases intron usage of chr22:50281499-50281566(-) |
| scz | Brain_Cerebellum | rs111925803 | C | -0.5527894496917725 | risk allele C decreases intron usage of chr22:50281499-50281643(-) |
| scz | Brain_Cerebellum | rs111925803 | C | -0.36653754115104675 | risk allele C decreases intron usage of chr22:50290597-50307553(-) |
| scz | Brain_Cerebellum | rs111925803 | C | 0.9847449660301208 | risk allele C increases intron usage of chr22:50294778-50295377(-) |
| scz | Brain_Cerebellum | rs111925803 | C | -0.663628876209259 | risk allele C decreases intron usage of chr22:50294778-50301347(-) |
| scz | Brain_Cerebellum | rs111925803 | C | -0.829789400100708 | risk allele C decreases intron usage of chr22:50294778-50307553(-) |
| scz | Brain_Cerebellum | rs111925803 | C | 0.6122559905052185 | risk allele C increases intron usage of chr22:50295429-50301347(-) |
| scz | Brain_Cerebellum | rs111925803 | C | 0.8459414839744568 | risk allele C increases intron usage of chr22:50295429-50307553(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, G3BP2, RBM14, CELF5, PABPC3, IGF2BP2
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ACO1(1), AGO1(1), AGO2(1), ANKHD1(1), CELF2(1), CELF4(1), CELF5(1), CELF6(1), CNOT4(1), CPEB2(1), CPEB4(1), CSTF2(1), DHX58(1), ELAVL1(1), ELAVL4(1), FXR1(1), FXR2(1), G3BP2(1), HNRNPA3(1)

## 5. Interpretation
PLXNB2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
