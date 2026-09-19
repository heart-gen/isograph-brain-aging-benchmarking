# TRAF3 — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** SCZ  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.312  ·  missense o/e 0.57

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.02 | no | — | yes | — | risk allele A increases TRAF3 expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Spinal_cord_cervical_c-1 | rs12588339 | A | 0.3340833783149719 | risk allele A increases expression of TRAF3 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SNRPB2, G3BP2, IGF2BP1, CNOT4, CELF4, CELF5, SNRNP70, PABPC3, YTHDC1, EIF4B, RBMS1, CPEB2, SAMD4A, ZNF638, IGF2BP2
- **All switched-motif RBPs (recurrence across regions):** CNOT4(7), SUPV3L1(7), PUM2(7), SNRNP70(7), ZNF638(7), EIF4B(7), CPEB2(1), CELF5(1), CELF4(1), ANKHD1(1), AKAP1(1), CSTF2(1), ESRP2(1), FXR2(1), IGF2BP1(1), IGF2BP2(1), KHDRBS2(1), MSI1(1), G3BP2(1), GRSF1(1)

## 5. Interpretation
TRAF3 colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
