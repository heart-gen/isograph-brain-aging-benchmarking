# PRSS36 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.13 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 1.218  ·  missense o/e 1.10

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.13 | no | — | no | — | risk allele G decreases PRSS36 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Nucleus_accumbens_basal_ganglia | rs78924645 | G | -0.7745634317398071 | risk allele G decreases expression of PRSS36 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** A1CF(1), ADAR(1), AGO1(1), CPEB2(1), CPEB4(1), CSTF2(1), DAZAP1(1), DDX58(1), DHX58(1), DHX9(1), ELAVL1(1), ELAVL2(1), ELAVL4(1), ENOX1(1), ESRP1(1), F2(1), HNRNPA0(1), HNRNPC(1), HNRNPD(1), HNRNPDL(1)

## 5. Interpretation
PRSS36 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
