# RANBP10 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** ALS  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.482  ·  missense o/e 0.82

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Cerebellum | 0.03 | no | — | no | — | risk allele G increases RANBP10 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellum | rs7195415 | G | 0.8765493035316467 | risk allele G increases expression of RANBP10 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** ACO1(1), AGO2(1), AKAP1(1), CELF4(1), CELF5(1), CNOT4(1), CPEB1(1), CPEB2(1), CPEB4(1), CSTF2(1), DDX19B(1), EIF4A3(1), ENOX1(1), ESRP1(1), G3BP1(1), GRSF1(1), HNRNPCL1(1), HNRNPDL(1), HNRNPLL(1), HNRNPU(1)

## 5. Interpretation
RANBP10 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
