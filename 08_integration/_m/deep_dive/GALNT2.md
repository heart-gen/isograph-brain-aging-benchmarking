# GALNT2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.07 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.515  ·  missense o/e 0.81

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cerebellar_Hemisphere | 0.07 | no | — | yes | — | risk allele A increases GALNT2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellar_Hemisphere | rs11807834 | A | 0.5212666392326355 | risk allele A increases expression of GALNT2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM41, RBMS1, RBM14, PUM1, RBMS3, G3BP1, SNRPB2, IGF2BP1, FXR1, ENOX1, PABPC3, ZCRB1, PABPN1, CELF6, AKAP1
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), AGO1(5), AGO2(5), AKAP1(5), CELF4(5), CELF5(5), CELF6(5), CNOT4(5), CPEB1(5), CPEB2(5), CPEB4(5), CSTF2(5), DAZAP1(5), DDX19B(5), DHX58(5), EIF4A3(5), ELAVL1(5), ELAVL3(5), ELAVL4(5), ENOX1(5)

## 5. Interpretation
GALNT2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
