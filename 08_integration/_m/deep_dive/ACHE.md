# ACHE — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.486  ·  missense o/e 0.73

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele C decreases usage of junction chr7:100894252-100894536(-) (unmapped transcript |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele C increases usage of junction chr7:100894252-100896380(-) (unmapped transcript |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele C decreases usage of junction chr7:100894255-100894536(-) (unmapped transcript |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele C increases usage of junction chr7:100894255-100896380(-) (unmapped transcript |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | no | no | — | risk allele C decreases usage of junction chr7:100896281-100896553(-) (ENST00000441605.2;  |
| ad | sQTL | Brain_Cerebellum | 0.01 | no | — | no | — | risk allele C increases usage of junction chr7:100896452-100896553(-) (unmapped transcript |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellum | rs3808356 | C | -0.6588960886001587 | risk allele C decreases intron usage of chr7:100894252-100894536(-) |
| ad | Brain_Cerebellum | rs3808356 | C | 0.869911253452301 | risk allele C increases intron usage of chr7:100894252-100896380(-) |
| ad | Brain_Cerebellum | rs3808356 | C | -0.917330801486969 | risk allele C decreases intron usage of chr7:100894255-100894536(-) |
| ad | Brain_Cerebellum | rs3808356 | C | 0.6134864687919617 | risk allele C increases intron usage of chr7:100894255-100896380(-) |
| ad | Brain_Cerebellum | rs3808356 | C | -0.8317129611968994 | risk allele C decreases intron usage of chr7:100896281-100896553(-) |
| ad | Brain_Cerebellum | rs3808356 | C | 0.8650103211402893 | risk allele C increases intron usage of chr7:100896452-100896553(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, G3BP2, AKAP1, G3BP1
- **All switched-motif RBPs (recurrence across regions):** AGO1(2), AKAP1(2), CELF6(2), CPEB1(2), DDX19B(2), EIF4B(2), ELAVL3(2), ENOX1(2), ESRP1(2), G3BP1(2), G3BP2(2), GRSF1(2), HNRNPA0(2), HNRNPAB(2), HNRNPC(2), HNRNPD(2), HNRNPDL(2), HNRNPU(2), IGF2BP1(2), KHDRBS1(2)

## 5. Interpretation
ACHE has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
