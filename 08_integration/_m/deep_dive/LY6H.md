# LY6H — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.03 (Brain_Hypothalamus)
- **Constraint:** LOEUF 1.081  ·  missense o/e 0.96

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellum | 0.02 | no | — | no | — | risk allele A increases LY6H expression (gene-level; no intron) |
| scz | sQTL | Brain_Hypothalamus | 0.03 | no | — | no | — | risk allele T increases usage of junction chr8:143159709-143160546(-) (unmapped transcript |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs7830479 | A | 0.6623439788818359 | risk allele A increases expression of LY6H |
| scz | Brain_Hypothalamus | rs11787024 | T | 0.8551557660102844 | risk allele T increases intron usage of chr8:143159709-143160546(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** CELF5, PTBP2
- **All switched-motif RBPs (recurrence across regions):** CELF4(3), CELF5(3), CPEB2(3), CSTF2(3), DHX58(3), ESRP1(3), GRSF1(3), HNRNPA0(3), HNRNPA2B1(3), HNRNPDL(3), NONO(3), NUDT21(3), PABPN1(3), PTBP2(3), RBFOX2(3), RBM24(3), RBM6(3), TRA2A(3)

## 5. Interpretation
LY6H has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
