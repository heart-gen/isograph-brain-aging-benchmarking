# FNBP1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.456  ·  missense o/e 0.82

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.02 | no | no | yes | — | risk allele G increases usage of junction chr9:129929695-129957360(-) (ENST00000355681.3,E |
| als | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.02 | no | no | yes | — | risk allele G decreases usage of junction chr9:129957464-129958491(-) (ENST00000355681.3,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Nucleus_accumbens_basal_ganglia | rs17519205 | G | 1.614669680595398 | risk allele G increases intron usage of chr9:129929695-129957360(-) |
| als | Brain_Nucleus_accumbens_basal_ganglia | rs17519205 | G | -1.6855639219284058 | risk allele G decreases intron usage of chr9:129957464-129958491(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS3, PPRC1, IGF2BP1, ZC3H10, HNRNPA3, CELF4, CELF5, CELF6, PABPC3, RBM46, PABPC5, YTHDC1, A1CF, ACO1, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), CELF4(5), CELF6(5), CELF5(5), CPEB2(5), CPEB4(5), EIF4A3(5), DDX19B(5), ESRP2(5), FXR2(5), PUM1(5), PPRC1(5), IGF2BP1(5), KHDRBS2(5), HNRNPCL1(5), HNRNPA3(5), SUPV3L1(5), U2AF2(5), YTHDC1(5)

## 5. Interpretation
FNBP1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
