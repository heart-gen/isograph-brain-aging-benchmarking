# EBPL — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level)

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.01 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 1.412  ·  missense o/e 0.99

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.01 | no | — | no | — | risk allele C decreases EBPL expression (gene-level; no intron) |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Spinal_cord_cervical_c-1 | rs117523952 | C | -1.2308763265609741 | risk allele C decreases expression of EBPL |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** HNRNPD, HNRNPA0, PABPC1, ZFP36, QKI
- **All switched-motif RBPs (recurrence across regions):** ACO1(1), AGO2(1), GRSF1(1), HNRNPA0(1), HNRNPAB(1), HNRNPCL1(1), HNRNPD(1), HNRNPDL(1), HNRNPM(1), IFIH1(1), IGHMBP2(1), NONO(1), OAS1(1), PABPC1(1), PABPC4(1), PUM1(1), QKI(1), RALY(1), RBM4(1), RBM41(1)

## 5. Interpretation
EBPL colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
