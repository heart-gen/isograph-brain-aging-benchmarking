# GSTO1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.12 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 1.236  ·  missense o/e 0.78

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.12 | no | — | yes | — | risk allele A decreases GSTO1 expression (gene-level; no intron) |
| scz | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.12 | no | — | no | — | risk allele A decreases GSTO1 expression (gene-level; no intron) |
| scz | sQTL | Brain_Cerebellum | 0.07 | no | no | yes | — | risk allele G decreases usage of junction chr10:104259798-104266084(+) (ENST00000369710.8, |
| scz | sQTL | Brain_Cerebellum | 0.07 | no | no | no | — | risk allele G decreases usage of junction chr10:104259798-104266084(+) (ENST00000369710.8, |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs10883990 | A | -0.22588178515434265 | risk allele A decreases expression of GSTO1 |
| scz | Brain_Cerebellum | rs17883150 | G | -0.5410272479057312 | risk allele G decreases intron usage of chr10:104259798-104266084(+) |
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs10883990 | A | -0.22588178515434265 | risk allele A decreases expression of GSTO1 |
| scz | Brain_Cerebellum | rs17883150 | G | -0.5410272479057312 | risk allele G decreases intron usage of chr10:104259798-104266084(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** AGO1(3), CPEB1(3), CSTF2(3), ELAVL3(3), ELAVL4(3), G3BP2(3), GRSF1(3), HNRNPA0(3), HNRNPDL(3), HNRNPM(3), IGF2BP3(3), NELFE(3), PABPC3(3), PHAX(3), PABPN1(3), PUM1(3), PTBP2(3), SRSF4(3), SRSF5(3), PUM2(3)

## 5. Interpretation
GSTO1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
