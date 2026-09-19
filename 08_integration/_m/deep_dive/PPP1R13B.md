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
- **Switched *and* module-enriched (q<0.05) RBPs:** DDX58, PABPC5, ZNF346, ADAR, IFIH1, RBM41, ERI1, PABPC4, CNOT4, MBNL1, RBMY1A1, CPEB2
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), ACO1(3), AGO1(3), AGO2(3), CELF5(3), CELF4(3), ELAVL3(3), EIF4B(3), CELF6(3), CNOT4(3), CPEB2(3), CSTF2(3), DHX58(3), DDX19B(3), EIF4A3(3), HNRNPM(3), HNRNPA3(3), HNRNPAB(3), HNRNPA0(3), G3BP1(3)

## 5. Interpretation
PPP1R13B colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
