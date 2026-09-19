# TMEM163 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.08 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 1.024  ·  missense o/e 0.85

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.08 | no | — | no | — | risk allele T increases TMEM163 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Nucleus_accumbens_basal_ganglia | rs6741007 | T | 0.3094671964645386 | risk allele T increases expression of TMEM163 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** A1CF, RBMS3, AGO1, CELF5, PABPC3, AKAP1, YTHDC1, PTBP2, TARDBP
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), AGO1(4), AKAP1(4), ANKHD1(4), CELF4(4), CELF5(4), CELF6(4), CSTF2(4), ESRP1(4), ESRP2(4), HNRNPA1L2(4), HNRNPA3(4), HNRNPM(4), IGF2BP1(4), KHDRBS2(4), KHDRBS3(4), LIN28A(4), NELFE(4), PABPC3(4)

## 5. Interpretation
TMEM163 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
