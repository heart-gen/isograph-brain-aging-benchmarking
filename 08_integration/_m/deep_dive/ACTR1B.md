# ACTR1B — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.56 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.968  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellum | 0.19 | no | — | yes | — | risk allele G decreases ACTR1B expression (gene-level; no intron) |
| scz | eQTL | Brain_Cerebellum | 0.18 | no | — | no | — | risk allele G decreases ACTR1B expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellum | 0.56 | no | no | yes | — | risk allele G decreases usage of junction chr2:97658316-97658427(-) (ENST00000289228.7,ENS |
| scz | sQTL | Brain_Cerebellum | 0.56 | no | — | yes | — | risk allele G increases usage of junction chr2:97658316-97658879(-) (unmapped transcript;  |
| scz | sQTL | Brain_Cerebellum | 0.56 | no | no | yes | — | risk allele G increases usage of junction chr2:97658643-97658879(-) (ENST00000289228.7; no |
| scz | sQTL | Brain_Cerebellum | 0.53 | no | no | no | — | risk allele G decreases usage of junction chr2:97658316-97658427(-) (ENST00000289228.7,ENS |
| scz | sQTL | Brain_Cerebellum | 0.53 | no | — | no | — | risk allele G increases usage of junction chr2:97658316-97658879(-) (unmapped transcript;  |
| scz | sQTL | Brain_Cerebellum | 0.53 | no | no | no | — | risk allele G increases usage of junction chr2:97658643-97658879(-) (ENST00000289228.7; no |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs11692435 | G | -1.1622798442840576 | risk allele G decreases expression of ACTR1B |
| scz | Brain_Cerebellum | rs11692435 | G | -1.7143243551254272 | risk allele G decreases intron usage of chr2:97658316-97658427(-) |
| scz | Brain_Cerebellum | rs11692435 | G | 1.0806782245635986 | risk allele G increases intron usage of chr2:97658316-97658879(-) |
| scz | Brain_Cerebellum | rs11692435 | G | 1.5108572244644165 | risk allele G increases intron usage of chr2:97658643-97658879(-) |
| scz | Brain_Cerebellum | rs11692435 | G | -1.1622798442840576 | risk allele G decreases expression of ACTR1B |
| scz | Brain_Cerebellum | rs11692435 | G | -1.7143243551254272 | risk allele G decreases intron usage of chr2:97658316-97658427(-) |
| scz | Brain_Cerebellum | rs11692435 | G | 1.0806782245635986 | risk allele G increases intron usage of chr2:97658316-97658879(-) |
| scz | Brain_Cerebellum | rs11692435 | G | 1.5108572244644165 | risk allele G increases intron usage of chr2:97658643-97658879(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, PABPC4, G3BP2, HNRNPLL, IGF2BP2, YTHDC1, RBM46
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), AGO1(5), CELF6(5), CPEB1(5), CPEB2(5), CPEB4(5), DHX58(5), ELAVL1(5), ELAVL3(5), ELAVL4(5), G3BP2(5), GRSF1(5), HNRNPA0(5), HNRNPA1L2(5), HNRNPA3(5), HNRNPC(5), HNRNPD(5), HNRNPLL(5), HNRNPU(5)

## 5. Interpretation
ACTR1B has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
