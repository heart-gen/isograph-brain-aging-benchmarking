# SEC61A2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.657  ·  missense o/e 0.55

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele T decreases usage of junction chr10:12158105-12160930(+) (ENST00000298428.14,E |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T increases usage of junction chr10:12158499-12160930(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | — | yes | — | risk allele T increases usage of junction chr10:12159626-12160930(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele T increases usage of junction chr10:12162289-12164268(+) (ENST00000298428.14,E |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele T decreases usage of junction chr10:12162289-12167753(+) (ENST00000475268.5,EN |
| scz | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele T decreases usage of junction chr10:12167893-12169258(+) (ENST00000475268.5,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs11257565 | T | -1.1014289855957031 | risk allele T decreases intron usage of chr10:12158105-12160930(+) |
| scz | Brain_Cerebellum | rs11257565 | T | 1.0239753723144531 | risk allele T increases intron usage of chr10:12158499-12160930(+) |
| scz | Brain_Cerebellum | rs11257565 | T | 0.8371362686157227 | risk allele T increases intron usage of chr10:12159626-12160930(+) |
| scz | Brain_Cerebellum | rs11257565 | T | 0.6068113446235657 | risk allele T increases intron usage of chr10:12162289-12164268(+) |
| scz | Brain_Cerebellum | rs11257565 | T | -0.47997769713401794 | risk allele T decreases intron usage of chr10:12162289-12167753(+) |
| scz | Brain_Cerebellum | rs11257565 | T | -0.7785850763320923 | risk allele T decreases intron usage of chr10:12167893-12169258(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO1(3), AGO2(3), CPEB2(3), CSTF2(3), DDX19B(3), DHX58(3), EIF4A3(3), EIF4B(3), ESRP1(3), FXR1(3), G3BP1(3), HNRNPA3(3), HNRNPLL(3), HNRNPM(3), HNRNPU(3), IGF2BP2(3), IGHMBP2(3), KHDRBS3(3), NELFE(3)

## 5. Interpretation
SEC61A2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
