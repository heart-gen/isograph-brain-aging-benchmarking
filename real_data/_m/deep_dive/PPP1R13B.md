# PPP1R13B — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.457  ·  missense o/e 0.84

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Caudate_basal_ganglia | 0.02 | no | — | yes | — | risk allele A increases PPP1R13B expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs66676135 | A | 0.18532906472682953 | risk allele A increases expression of PPP1R13B |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** KHDRBS1, HNRNPD, ELAVL3, RBM41, PPIE, G3BP1, CPEB4, TIAL1, CPEB1, HNRNPDL, RBM6, AGO1, DHX58, ZC3H10
- **All switched-motif RBPs (recurrence across regions):** A1CF(5), ACO1(5), AGO1(5), AGO2(5), CELF4(5), ELAVL3(5), CELF5(5), CELF6(5), CNOT4(5), CPEB2(5), CSTF2(5), DHX58(5), DDX19B(5), EIF4A3(5), EIF4B(5), HNRNPU(5), HNRNPA3(5), G3BP1(5), HNRNPA0(5), FXR1(5)

## 5. Interpretation
PPP1R13B colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
