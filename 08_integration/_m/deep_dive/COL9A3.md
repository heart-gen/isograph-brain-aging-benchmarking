# COL9A3 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.026  ·  missense o/e 1.15

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele T increases usage of junction chr20:62818553-62819222(+) (ENST00000452372.2,EN |
| ad | sQTL | Brain_Caudate_basal_ganglia | 0.01 | no | no | yes | — | risk allele T decreases usage of junction chr20:62819293-62819929(+) (ENST00000452372.2,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Caudate_basal_ganglia | rs61734651 | T | 1.1272915601730347 | risk allele T increases intron usage of chr20:62818553-62819222(+) |
| ad | Brain_Caudate_basal_ganglia | rs61734651 | T | -1.0607349872589111 | risk allele T decreases intron usage of chr20:62819293-62819929(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPA0, DAZAP1, HNRNPU, PIWIL1, HNRNPK, SYNCRIP, PABPC1, ZFP36, NUDT21, QKI, U2AF2, RBM4, ELAVL3, RBM5
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), AGO1(3), AGO2(3), AKAP1(3), CELF1(3), CELF2(3), CELF4(3), CELF5(3), CELF6(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DAZAP1(3), DDX19B(3), DHX58(3), EIF4A3(3), ENOX1(3), ESRP2(3), ESRP1(3)

## 5. Interpretation
COL9A3 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
