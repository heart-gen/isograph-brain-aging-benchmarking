# MARK2 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Cortex)
- **Constraint:** LOEUF 0.201  ·  missense o/e 0.57

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Cortex | 0.03 | no | — | yes | — | risk allele T increases MARK2 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cortex | rs7121067 | T | 0.1440276950597763 | risk allele T increases expression of MARK2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBMS3, IGF2BP1, HNRNPCL1, PABPC5, YBX2, ESRP1, G3BP2, ZFP36L2, IGF2BP2, RALY, PABPN1, PABPC3, PABPC4, ACO1, KHDRBS1
- **All switched-motif RBPs (recurrence across regions):** ACO1(5), CELF5(5), CPEB1(5), CNOT4(5), DDX19B(5), DHX58(5), FXR2(5), EIF4A3(5), HNRNPLL(5), HNRNPCL1(5), YBX2(5), U2AF2(5), PABPC4(5), KHDRBS3(5), KHDRBS2(5), KHDRBS1(5), ZFP36L2(5), RALY(5), PABPN1(4), PABPC5(4)

## 5. Interpretation
MARK2 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
