# NUP85 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** ALS  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.329  ·  missense o/e 0.89

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A increases usage of junction chr17:75215823-75216290(+) (ENST00000577208.5; n |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A increases usage of junction chr17:75215823-75216311(+) (ENST00000583548.1; n |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A decreases usage of junction chr17:75215823-75218185(+) (ENST00000245544.9,EN |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele A increases usage of junction chr17:75216500-75218185(+) (unmapped transcript; |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele A increases usage of junction chr17:75216500-75225103(+) (unmapped transcript; |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A decreases usage of junction chr17:75226157-75226745(+) (ENST00000449421.6; n |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A increases usage of junction chr17:75226157-75228998(+) (ENST00000580879.5; n |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele A increases usage of junction chr17:75226157-75230953(+) (unmapped transcript; |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A increases usage of junction chr17:75226157-75231340(+) (ENST00000245544.9,EN |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A decreases usage of junction chr17:75228384-75228998(+) (ENST00000581335.1; n |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele A decreases usage of junction chr17:75229131-75230546(+) (unmapped transcript; |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele A decreases usage of junction chr17:75229131-75230953(+) (unmapped transcript; |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | no | yes | — | risk allele A decreases usage of junction chr17:75229131-75231340(+) (ENST00000540768.5,EN |
| als | sQTL | Brain_Cerebellar_Hemisphere | 0.02 | no | — | yes | — | risk allele A decreases usage of junction chr17:75230636-75231340(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | 1.0538437366485596 | risk allele A increases intron usage of chr17:75215823-75216290(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | 0.6520705223083496 | risk allele A increases intron usage of chr17:75215823-75216311(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | -0.7405231595039368 | risk allele A decreases intron usage of chr17:75215823-75218185(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | 0.5394098162651062 | risk allele A increases intron usage of chr17:75216500-75218185(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | 0.34252437949180603 | risk allele A increases intron usage of chr17:75216500-75225103(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | -1.0492968559265137 | risk allele A decreases intron usage of chr17:75226157-75226745(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | 1.2281837463378906 | risk allele A increases intron usage of chr17:75226157-75228998(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | 0.5850234627723694 | risk allele A increases intron usage of chr17:75226157-75230953(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | 1.2424386739730835 | risk allele A increases intron usage of chr17:75226157-75231340(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | -0.9058374166488647 | risk allele A decreases intron usage of chr17:75228384-75228998(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | -0.7091026902198792 | risk allele A decreases intron usage of chr17:75229131-75230546(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | -0.647588312625885 | risk allele A decreases intron usage of chr17:75229131-75230953(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | -1.114695429801941 | risk allele A decreases intron usage of chr17:75229131-75231340(+) |
| als | Brain_Cerebellar_Hemisphere | rs1478785 | A | -0.6740987300872803 | risk allele A decreases intron usage of chr17:75230636-75231340(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, G3BP2, AKAP1, SNRPB2
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), ADAR(3), AGO2(3), AKAP1(3), CNOT4(3), CPEB2(3), CSTF2(3), EIF4A3(3), DDX19B(3), DHX58(3), HNRNPCL1(3), HNRNPAB(3), HNRNPA3(3), HNRNPA0(3), G3BP2(3), LIN28A(3), KHDRBS3(3), IGF2BP2(3), KHDRBS2(3), HNRNPU(3)

## 5. Interpretation
NUP85 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
