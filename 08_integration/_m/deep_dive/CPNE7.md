# CPNE7 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Putamen_basal_ganglia)
- **Constraint:** LOEUF 1.522  ·  missense o/e 1.41

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Putamen_basal_ganglia | 0.01 | no | — | no | — | risk allele A increases CPNE7 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Putamen_basal_ganglia | rs4785751 | A | 0.1825851947069168 | risk allele A increases expression of CPNE7 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPLL, RBM14, IGF2BP2, PTBP2, ADAR
- **All switched-motif RBPs (recurrence across regions):** ACO1(3), ADAR(3), AGO1(3), AGO2(3), CELF2(3), CPEB1(3), CPEB4(3), CPEB2(3), CSTF2(3), DDX19B(3), ESRP1(3), DHX58(3), EIF4A3(3), ELAVL1(3), ELAVL2(3), ELAVL3(3), ELAVL4(3), ENOX1(3), HNRNPAB(3), FXR1(3)

## 5. Interpretation
CPNE7 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
