# NPEPL1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 1.225  ·  missense o/e 1.06

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele T decreases usage of junction chr20:58701158-58705498(+) (ENST00000533788.1; n |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele T increases usage of junction chr20:58701158-58707123(+) (ENST00000356091.11,E |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele T decreases usage of junction chr20:58705403-58705498(+) (unmapped transcript; |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele T decreases usage of junction chr20:58705569-58707123(+) (ENST00000533788.1; n |
| scz | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele T decreases usage of junction chr20:58713543-58713961(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs35104993 | T | -0.9315313696861267 | risk allele T decreases intron usage of chr20:58701158-58705498(+) |
| scz | Brain_Cerebellum | rs35104993 | T | 0.9402583241462708 | risk allele T increases intron usage of chr20:58701158-58707123(+) |
| scz | Brain_Cerebellum | rs35104993 | T | -0.694797933101654 | risk allele T decreases intron usage of chr20:58705403-58705498(+) |
| scz | Brain_Cerebellum | rs35104993 | T | -0.9620973467826843 | risk allele T decreases intron usage of chr20:58705569-58707123(+) |
| scz | Brain_Cerebellum | rs35104993 | T | -1.336612582206726 | risk allele T decreases intron usage of chr20:58713543-58713961(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, G3BP1, SNRPB2, ENOX1, ZC3H10, IGF2BP1, HNRNPA3, CELF4, CELF5, CELF6, MATR3, SNRNP70, RBM46, PABPC5, YTHDC1
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), AGO1(2), CPEB4(2), CPEB2(2), ESRP1(2), EIF4B(2), EIF4A3(2), DHX58(2), DDX19B(2), ENOX1(2), HNRNPU(2), HNRNPCL1(2), HNRNPA3(2), HNRNPA0(2), G3BP1(2), FXR2(2), ESRP2(2), KHDRBS3(2), KHDRBS2(2), IGF2BP2(2)

## 5. Interpretation
NPEPL1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
