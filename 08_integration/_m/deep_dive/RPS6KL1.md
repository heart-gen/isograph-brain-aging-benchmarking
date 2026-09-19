# RPS6KL1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.09 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.978  ·  missense o/e 0.97

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele C decreases RPS6KL1 expression (gene-level; no intron) |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.09 | no | no | yes | — | risk allele A increases usage of junction chr14:74911841-74918513(-) (ENST00000354625.6,EN |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.09 | no | no | yes | — | risk allele A decreases usage of junction chr14:74911841-74919845(-) (ENST00000555009.5; n |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.09 | no | no | yes | — | risk allele A increases usage of junction chr14:74918605-74919845(-) (ENST00000354625.6,EN |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.09 | no | no | yes | — | risk allele A increases usage of junction chr14:74921561-74922341(-) (ENST00000555647.5; n |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.09 | no | no | yes | — | risk allele A increases usage of junction chr14:74921561-74923180(-) (ENST00000555834.5; n |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs929579 | C | -0.15605716407299042 | risk allele C decreases expression of RPS6KL1 |
| als | Brain_Cerebellar_Hemisphere | rs7158047 | A | 0.29705825448036194 | risk allele A increases intron usage of chr14:74911841-74918513(-) |
| als | Brain_Cerebellar_Hemisphere | rs7158047 | A | -0.2730770409107208 | risk allele A decreases intron usage of chr14:74911841-74919845(-) |
| als | Brain_Cerebellar_Hemisphere | rs7158047 | A | 0.31739869713783264 | risk allele A increases intron usage of chr14:74918605-74919845(-) |
| als | Brain_Cerebellar_Hemisphere | rs7158047 | A | 0.32527413964271545 | risk allele A increases intron usage of chr14:74921561-74922341(-) |
| als | Brain_Cerebellar_Hemisphere | rs7158047 | A | 0.6292232871055603 | risk allele A increases intron usage of chr14:74921561-74923180(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** IGF2BP3, KHDRBS1, TIAL1, SYNCRIP, PABPC1, ELAVL3
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO2(3), CELF4(3), CELF5(3), CELF6(3), CNOT4(3), CPEB1(3), CPEB2(3), CPEB4(3), CSTF2(3), DHX58(3), EIF4B(3), ELAVL3(3), ENOX1(3), FXR1(3), FXR2(3), G3BP1(3), HNRNPCL1(3), HNRNPLL(3)

## 5. Interpretation
RPS6KL1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
