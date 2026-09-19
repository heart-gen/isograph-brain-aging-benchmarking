# ALG12 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.965  ·  missense o/e 0.90

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Caudate_basal_ganglia | 0.01 | no | — | no | — | risk allele T increases ALG12 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs7287579 | T | 0.47082826495170593 | risk allele T increases expression of ALG12 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PUM1, SF1, ZRANB2, SRSF6, TARDBP, CPEB4
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), AGO2(2), CELF4(2), CELF5(2), CPEB1(2), CPEB4(2), CSTF2(2), DAZAP1(2), DDX19B(2), EIF4B(2), ELAVL2(2), ELAVL3(2), ESRP1(2), GRSF1(2), HNRNPA0(2), HNRNPAB(2), HNRNPC(2), HNRNPCL1(2), HNRNPD(2), HNRNPDL(2)

## 5. Interpretation
ALG12 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
