# TMEM219 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 1.425  ·  missense o/e 0.94

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.02 | no | — | yes | — | risk allele G increases TMEM219 expression (gene-level; no intron) |
| scz | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.02 | no | — | no | — | risk allele G increases TMEM219 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs9932196 | G | 0.09732464700937271 | risk allele G increases expression of TMEM219 |
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs9932196 | G | 0.09732464700937271 | risk allele G increases expression of TMEM219 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** ZFP36L2, PABPC5, IFIH1, CNOT4, RBM46
- **All switched-motif RBPs (recurrence across regions):** ACO1(2), AGO2(2), CNOT4(2), CSTF2(2), DAZAP1(2), EIF4B(2), ESRP1(2), HNRNPCL1(2), HNRNPA0(2), SRSF7(2), SUPV3L1(2), HNRNPD(2), HNRNPDL(2), IGHMBP2(2), IFIH1(2), KHDRBS2(2), KHDRBS3(2), PCBP1(2), LIN28A(2), PABPC1(2)

## 5. Interpretation
TMEM219 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
