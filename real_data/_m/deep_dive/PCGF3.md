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
- **Switched *and* module-enriched (q<0.05) RBPs:** SART3, YTHDC1, SUPV3L1, SNRPB2, AGO1, ZFP36L2, AKAP1, FXR1, AGO2, RBMS3, HNRNPLL, CELF4, A1CF, ENOX1, CELF5
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), AGO1(5), AGO2(5), AKAP1(5), CELF4(5), CELF5(5), CPEB2(5), CSTF2(5), DHX58(5), ENOX1(5), ESRP1(5), FXR1(5), FXR2(5), G3BP2(5), HNRNPA0(5), HNRNPA1L2(5), HNRNPA3(5), HNRNPCL1(5), HNRNPLL(5)

## 5. Interpretation
PCGF3 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
