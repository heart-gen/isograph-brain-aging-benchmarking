# DUS2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 0.859  ·  missense o/e 0.86

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Putamen_basal_ganglia | 0.02 | no | — | no | — | risk allele C decreases DUS2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Putamen_basal_ganglia | rs60342335 | C | -0.7407863736152649 | risk allele C decreases expression of DUS2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMY1A1
- **All switched-motif RBPs (recurrence across regions):** ACO1(1), ADAR(1), AGO1(1), CPEB1(1), CPEB4(1), CSTF2(1), DAZAP1(1), DDX19B(1), DHX9(1), EIF4A3(1), EIF4B(1), ELAVL3(1), ENOX1(1), ESRP2(1), FXR1(1), GRSF1(1), HNRNPA0(1), HNRNPA2B1(1), HNRNPA3(1), HNRNPAB(1)

## 5. Interpretation
DUS2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
