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
- **Switched *and* module-enriched (q<0.05) RBPs:** RC3H1, PPIE, OAS1, HNRNPA0, RBM14, AKAP1, FXR1, ZFP36L2, SUPV3L1, CELF4, RBM24, IGF2BP3
- **All switched-motif RBPs (recurrence across regions):** ACO1(4), AGO1(4), AGO2(4), AKAP1(4), CELF4(4), CELF5(4), CNOT4(4), CPEB2(4), EIF4A3(4), ENOX1(4), ESRP1(4), FXR1(4), G3BP1(4), G3BP2(4), HNRNPA0(4), HNRNPA1L2(4), HNRNPA3(4), HNRNPAB(4), HNRNPDL(4), HNRNPLL(4)

## 5. Interpretation
PAK6 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
