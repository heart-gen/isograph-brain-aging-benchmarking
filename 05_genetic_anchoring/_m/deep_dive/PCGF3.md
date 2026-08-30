# PCGF3 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** AD,PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.445  ·  missense o/e 0.51

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | no | yes | — | risk allele A increases usage of junction chr4:705970-730630(+) (ENST00000362003.10,ENST00 |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | no | yes | — | risk allele A decreases usage of junction chr4:725229-730630(+) (ENST00000433814.5; not in |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | — | yes | — | risk allele A decreases usage of junction chr4:728628-730630(+) (unmapped transcript; not  |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | no | yes | — | risk allele A increases usage of junction chr4:731110-733672(+) (ENST00000362003.10,ENST00 |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | no | yes | — | risk allele A decreases usage of junction chr4:732498-733672(+) (ENST00000419774.5,ENST000 |
| pd | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr4:731110-733672(+) (ENST00000362003.10,ENST00 |
| pd | sQTL | Brain_Cerebellum | 0.01 | no | no | yes | — | risk allele G increases usage of junction chr4:761416-764984(+) (ENST00000362003.10,ENST00 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PUM1
- **All switched-motif RBPs (recurrence across regions):** ACO1(6), AGO1(6), AGO2(6), AKAP1(6), CELF4(6), CELF5(6), CELF6(6), CNOT4(6), CPEB2(6), CSTF2(6), DHX58(6), ENOX1(6), ESRP1(6), FXR1(6), FXR2(6), G3BP2(6), HNRNPA0(6), HNRNPA1L2(6), HNRNPA3(6), HNRNPLL(6)

## 5. Interpretation
PCGF3 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
