# ZSWIM7 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.02 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 1.603  ·  missense o/e 1.16

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C decreases usage of junction chr17:15977092-15977579(-) (ENST00000490395.5,EN |
| pd | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C decreases usage of junction chr17:15977092-15978044(-) (ENST00000476496.5; n |
| pd | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C increases usage of junction chr17:15977913-15978044(-) (ENST00000460315.5,EN |
| pd | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C decreases usage of junction chr17:15978163-15980243(-) (ENST00000497434.6; n |
| pd | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C increases usage of junction chr17:15993778-15999519(-) (ENST00000399277.6,EN |
| pd | sQTL | Brain_Caudate_basal_ganglia | 0.02 | no | no | yes | — | risk allele C decreases usage of junction chr17:15993778-15999660(-) (ENST00000399280.6,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Caudate_basal_ganglia | rs1045599 | C | -1.3808925151824951 | risk allele C decreases intron usage of chr17:15977092-15977579(-) |
| pd | Brain_Caudate_basal_ganglia | rs1045599 | C | -0.46644091606140137 | risk allele C decreases intron usage of chr17:15977092-15978044(-) |
| pd | Brain_Caudate_basal_ganglia | rs1045599 | C | 1.4173133373260498 | risk allele C increases intron usage of chr17:15977913-15978044(-) |
| pd | Brain_Caudate_basal_ganglia | rs1045599 | C | -0.5074920654296875 | risk allele C decreases intron usage of chr17:15978163-15980243(-) |
| pd | Brain_Caudate_basal_ganglia | rs1045599 | C | 0.38954389095306396 | risk allele C increases intron usage of chr17:15993778-15999519(-) |
| pd | Brain_Caudate_basal_ganglia | rs1045599 | C | -0.4231342077255249 | risk allele C decreases intron usage of chr17:15993778-15999660(-) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPC, HNRNPA0, DAZAP1, HNRNPU, HNRNPK, SYNCRIP, U2AF2, ELAVL3
- **All switched-motif RBPs (recurrence across regions):** AGO1(2), CNOT4(2), CPEB2(2), CPEB4(2), DAZAP1(2), DDX19B(2), DHX58(2), EIF4A3(2), EIF4B(2), ELAVL3(2), ENOX1(2), ESRP1(2), GRSF1(2), HNRNPA0(2), HNRNPA2B1(2), HNRNPA3(2), HNRNPC(2), HNRNPCL1(2), HNRNPD(2), HNRNPDL(2)

## 5. Interpretation
ZSWIM7 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
