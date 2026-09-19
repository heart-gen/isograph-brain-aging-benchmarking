# ITGB1BP1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.06 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.905  ·  missense o/e 0.74

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.06 | no | no | no | — | risk allele C increases usage of junction chr2:9412405-9418626(-) (ENST00000464228.5; not  |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.06 | no | — | no | — | risk allele C decreases usage of junction chr2:9418732-9418818(-) (unmapped transcript; no |
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.06 | no | no | no | — | risk allele C increases usage of junction chr2:9418732-9419990(-) (ENST00000360635.7,ENST0 |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellar_Hemisphere | rs10169262 | C | 0.37140506505966187 | risk allele C increases intron usage of chr2:9412405-9418626(-) |
| ad | Brain_Cerebellar_Hemisphere | rs10169262 | C | -0.4162002503871918 | risk allele C decreases intron usage of chr2:9418732-9418818(-) |
| ad | Brain_Cerebellar_Hemisphere | rs10169262 | C | 0.5649600625038147 | risk allele C increases intron usage of chr2:9418732-9419990(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPU, HNRNPK, TIAL1, SYNCRIP, U2AF2
- **All switched-motif RBPs (recurrence across regions):** ADAR(2), CELF4(2), AGO2(2), AKAP1(2), ANKHD1(2), CNOT4(2), CELF5(2), CPEB2(2), CPEB4(2), ZNF638(2), CSTF2(2), DDX19B(2), DHX58(2), DHX9(2), EIF4A3(2), ENOX1(2), ESRP1(2), ESRP2(2), G3BP1(2), HNRNPA3(2)

## 5. Interpretation
ITGB1BP1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
