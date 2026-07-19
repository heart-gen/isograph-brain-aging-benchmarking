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
- **Switched *and* module-enriched (q<0.05) RBPs:** SART3, A1CF, CPEB2, YBX2, YTHDC1, RALY, DHX58, SNRPB2, HNRNPAB, ZFP36L2, AKAP1, FXR1, AGO2, CELF4, PABPC4
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), AGO2(4), AKAP1(4), CELF4(4), CELF5(4), CELF6(4), CPEB1(4), CPEB2(4), DDX19B(4), DHX58(4), ESRP1(4), EIF4A3(4), EIF4B(4), ELAVL3(4), FXR2(4), FXR1(4), G3BP1(4), HNRNPA0(4), RBFOX2(4)

## 5. Interpretation
PPP1R13B colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
