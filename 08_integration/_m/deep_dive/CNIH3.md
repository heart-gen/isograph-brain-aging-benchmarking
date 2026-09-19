# CNIH3 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.984  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224434862-224454278(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T increases usage of junction chr1:224434862-224456854(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224439774-224454278(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224454343-224456854(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224454343-224456859(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224456921-224459124(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224456921-224459143(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224456921-224459210(+) (unmapped transcript |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | no | — | risk allele T decreases usage of junction chr1:224456921-224459213(+) (unmapped transcript |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -0.6716414093971252 | risk allele T decreases intron usage of chr1:224434862-224454278(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | 0.6768707036972046 | risk allele T increases intron usage of chr1:224434862-224456854(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -0.6070449352264404 | risk allele T decreases intron usage of chr1:224439774-224454278(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -0.5411155223846436 | risk allele T decreases intron usage of chr1:224454343-224456854(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -0.5272606015205383 | risk allele T decreases intron usage of chr1:224454343-224456859(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -0.4843730330467224 | risk allele T decreases intron usage of chr1:224456921-224459124(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -0.9727810025215149 | risk allele T decreases intron usage of chr1:224456921-224459143(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -0.7607853412628174 | risk allele T decreases intron usage of chr1:224456921-224459210(+) |
| als | Brain_Cerebellar_Hemisphere | rs10799573 | T | -1.118112564086914 | risk allele T decreases intron usage of chr1:224456921-224459213(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, HNRNPA0, HNRNPC, ELAVL3, KHDRBS1, HNRNPD, A1CF, CPEB1, SNRPB2, RNASEL, SRSF10, PABPC1, CPEB4, PABPC4, U2AF2
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), ACO1(6), AGO1(6), AGO2(6), CNOT4(6), CELF6(6), CPEB1(6), CPEB2(6), DDX19B(6), CPEB4(6), CSTF2(6), DAZAP1(6), DHX58(6), DDX58(6), ENOX1(6), EIF4A3(6), U2AF2(6), YTHDC1(6), G3BP2(6), FXR1(6)

## 5. Interpretation
CNIH3 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
