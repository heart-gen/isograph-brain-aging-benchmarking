# DGKQ — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** ALS,LBD,PD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 1.00 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 1.155  ·  missense o/e 0.99

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | eQTL | Brain_Cerebellar_Hemisphere | 1.00 | no | — | no | — | risk allele C decreases DGKQ expression (gene-level; no intron) |
| als | eQTL | Brain_Cerebellar_Hemisphere | 0.10 | no | — | no | — | risk allele C decreases DGKQ expression (gene-level; no intron) |
| lbd | eQTL | Brain_Cerebellar_Hemisphere | 0.03 | no | — | no | — | risk allele C decreases DGKQ expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| als | Brain_Cerebellar_Hemisphere | rs34311866 | C | -0.2708304524421692 | risk allele C decreases expression of DGKQ |
| lbd | Brain_Cerebellar_Hemisphere | rs34311866 | C | -0.2708304524421692 | risk allele C decreases expression of DGKQ |
| pd | Brain_Cerebellar_Hemisphere | rs34311866 | C | -0.2708304524421692 | risk allele C decreases expression of DGKQ |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM6, SNRPB2, RBMS3, PPRC1, ENOX1, IGF2BP1, CELF4, CELF5, FXR1, CELF6, EIF4B, RBMS1, A1CF, ACO1, ZFP36L2
- **All switched-motif RBPs (recurrence across regions):** A1CF(4), ACO1(4), AGO1(4), AGO2(4), AKAP1(4), CELF2(4), CELF4(4), CELF5(4), CELF6(4), CSTF2(4), DAZAP1(4), DDX19B(4), DHX58(4), EIF4B(4), ELAVL1(4), ELAVL2(4), ELAVL3(4), ELAVL4(4), ENOX1(4), ESRP1(4)

## 5. Interpretation
DGKQ colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
