# ITSN1 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.03 (Brain_Caudate_basal_ganglia)
- **Constraint:** LOEUF 0.301  ·  missense o/e 0.78

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Caudate_basal_ganglia | 0.03 | no | — | yes | — | risk allele G increases ITSN1 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Caudate_basal_ganglia | rs2251854 | G | 0.2255571335554123 | risk allele G increases expression of ITSN1 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPDL, HNRNPA0, ELAVL3, KHDRBS1, HNRNPD, A1CF, RBMS3, RBM41, CPEB1, SNRPB2, SRSF10, PABPC1, RBMY1A1, G3BP2, AGO1
- **All switched-motif RBPs (recurrence across regions):** A1CF(6), AGO2(6), AKAP1(6), CELF5(6), CELF4(6), CELF6(6), CPEB1(6), FXR2(6), CPEB2(6), CSTF2(6), DDX19B(6), EIF4A3(6), DHX58(6), ELAVL3(6), EIF4B(6), G3BP1(6), FXR1(6), ENOX1(6), ESRP2(6), HNRNPAB(6)

## 5. Interpretation
ITSN1 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
