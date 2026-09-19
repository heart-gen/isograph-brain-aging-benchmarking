# ADAM10 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.04 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.276  ·  missense o/e 0.65

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.04 | no | — | yes | — | risk allele T increases ADAM10 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Spinal_cord_cervical_c-1 | rs1427281 | T | 0.3229811489582062 | risk allele T increases expression of ADAM10 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** PPRC1, RBM41
- **All switched-motif RBPs (recurrence across regions):** A1CF(2), ACO1(2), ADAR(2), AGO2(2), AKAP1(2), CELF4(2), DHX9(2), CELF5(2), CELF6(2), CNOT4(2), CPEB2(2), CSTF2(2), DDX19B(2), DHX58(2), G3BP1(2), EIF4A3(2), ESRP1(2), ESRP2(2), G3BP2(2), FXR1(2)

## 5. Interpretation
ADAM10 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
