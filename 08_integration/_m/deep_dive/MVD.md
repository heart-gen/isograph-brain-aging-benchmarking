# MVD — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 1.602  ·  missense o/e 1.22

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.01 | no | — | yes | — | risk allele G decreases MVD expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Nucleus_accumbens_basal_ganglia | rs12929177 | G | -0.23933951556682587 | risk allele G decreases expression of MVD |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, HNRNPC, ELAVL3, HNRNPD, ELAVL4, IGF2BP3, CPEB1, SRSF10, PABPC1, CPEB4, ZFP36, TIAL1, U2AF2, PPIE, RBMY1A1
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), AGO1(4), CPEB1(4), CELF4(4), CELF6(4), GRSF1(4), HNRNPA3(4), DDX19B(4), DAZAP1(4), CPEB4(4), HNRNPAB(4), HNRNPA2B1(4), ELAVL3(4), ELAVL2(4), G3BP1(4), IGHMBP2(4), IGF2BP3(4), IGF2BP2(4), IFIH1(4), HNRNPU(4)

## 5. Interpretation
MVD colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
