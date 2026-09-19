# SLCO1A2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 1.214  ·  missense o/e 0.93

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr12:21306988-21318782(-) (ENST00000544020.5,EN |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T increases usage of junction chr12:21318923-21319372(-) (ENST00000445053.1,EN |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr12:21318923-21324594(-) (ENST00000480394.5; n |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T decreases usage of junction chr12:21318923-21334588(-) (ENST00000307378.10,E |
| als | sQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | no | no | — | risk allele T increases usage of junction chr12:21319579-21324594(-) (ENST00000463718.5; n |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Spinal_cord_cervical_c-1 | rs2306229 | T | -0.4848131537437439 | risk allele T decreases intron usage of chr12:21306988-21318782(-) |
| als | Brain_Spinal_cord_cervical_c-1 | rs2306229 | T | 0.6213617324829102 | risk allele T increases intron usage of chr12:21318923-21319372(-) |
| als | Brain_Spinal_cord_cervical_c-1 | rs2306229 | T | -0.5239601731300354 | risk allele T decreases intron usage of chr12:21318923-21324594(-) |
| als | Brain_Spinal_cord_cervical_c-1 | rs2306229 | T | -0.4562588930130005 | risk allele T decreases intron usage of chr12:21318923-21334588(-) |
| als | Brain_Spinal_cord_cervical_c-1 | rs2306229 | T | 0.6317012906074524 | risk allele T increases intron usage of chr12:21319579-21324594(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, CPEB1, PABPC1, CPEB4, DAZAP1, HNRNPU, U2AF2, HNRNPK, SYNCRIP, PPIE, AGO1, PUM1, SF1, ZRANB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), ADAR(2), AGO1(2), AGO2(2), AKAP1(2), CELF5(2), CELF6(2), CNOT4(2), CPEB2(2), ESRP1(2), CPEB4(2), CSTF2(2), DAZAP1(2), DDX58(2), DDX19B(2), DHX9(2), EIF4A3(2), FUS(2), EIF4B(2)

## 5. Interpretation
SLCO1A2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
