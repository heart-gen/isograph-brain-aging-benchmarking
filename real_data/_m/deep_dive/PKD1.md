# PKD1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.360  ·  missense o/e 1.22

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele G decreases usage of junction chr16:2100566-2102102(-) (unmapped transcript; n |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele G decreases usage of junction chr16:2105662-2105865(-) (ENST00000564865.5; not |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele G increases usage of junction chr16:2105935-2106091(-) (unmapped transcript; n |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | no | yes | — | risk allele G decreases usage of junction chr16:2106024-2106091(-) (ENST00000262304.9,ENST |
| pd | sQTL | Brain_Cerebellum | 0.02 | no | — | yes | — | risk allele G decreases usage of junction chr16:2112399-2112788(-) (unmapped transcript; n |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** YBX2, CELF6, ZFP36L2, IGF2BP1, RBM24, FXR1, ENOX1, PABPN1, HNRNPA1L2, ESRP1, AKAP1, ESRP2, CNOT4, DHX58, CELF5
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), AGO1(4), AKAP1(4), CELF4(4), CELF5(4), CELF6(4), CMTR1(4), CNOT4(4), CPEB4(4), DDX19B(4), DHX58(4), ELAVL2(4), ELAVL3(4), ENOX1(4), ESRP1(4), ESRP2(4), FXR1(4), G3BP1(4), HNRNPA0(4), HNRNPA1L2(4)

## 5. Interpretation
PKD1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
