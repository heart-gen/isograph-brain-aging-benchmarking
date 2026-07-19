# PAK6 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.17 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.715  ·  missense o/e 0.92

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellum | 0.17 | no | — | no | — | risk allele C increases PAK6 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs56282503 | C | 0.35913631319999695 | risk allele C increases expression of PAK6 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM14, G3BP1, SAMD4A, RC3H1, RBM8A, OAS1, HNRNPA0, PPRC1, SUPV3L1, IGF2BP3, QKI
- **All switched-motif RBPs (recurrence across regions):** AGO1(3), AGO2(3), CELF4(3), CELF5(3), CPEB2(3), CPEB4(3), EIF4A3(3), ENOX1(3), ESRP1(3), G3BP1(3), G3BP2(3), HNRNPA0(3), HNRNPA1L2(3), HNRNPAB(3), HNRNPLL(3), IGF2BP3(3), NELFE(3), OAS1(3), PHAX(3), PPRC1(3)

## 5. Interpretation
PAK6 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
