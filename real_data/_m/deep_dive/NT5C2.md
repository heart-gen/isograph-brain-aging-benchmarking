# NT5C2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.128  ·  missense o/e 0.65

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.01 | no | — | yes | — | risk allele T decreases NT5C2 expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele G increases usage of junction chr10:103174982-103181185(-) (ENST00000404739.8, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs11191580 | T | -0.4148642420768738 | risk allele T decreases expression of NT5C2 |
| scz | Brain_Cerebellum | rs12412038 | G | 1.4276000261306763 | risk allele G increases intron usage of chr10:103174982-103181185(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ESRP1(3), SRSF11(3), PABPN1(3), PPRC1(3), RBM6(3), RBM8A(3), CELF4(1), DHX58(1), RBM41(1)

## 5. Interpretation
NT5C2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
