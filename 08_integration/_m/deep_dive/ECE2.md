# ECE2 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** PD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.10 (Brain_Nucleus_accumbens_basal_ganglia)
- **Constraint:** LOEUF 0.932  ·  missense o/e 0.88

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| pd | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.10 | no | no | no | — | risk allele C decreases usage of junction chr3:184276567-184276892(+) (ENST00000404464.8,E |
| pd | sQTL | Brain_Nucleus_accumbens_basal_ganglia | 0.10 | no | no | no | — | risk allele C increases usage of junction chr3:184276705-184276892(+) (ENST00000357474.9,E |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| pd | Brain_Nucleus_accumbens_basal_ganglia | rs843351 | C | -0.714849591255188 | risk allele C decreases intron usage of chr3:184276567-184276892(+) |
| pd | Brain_Nucleus_accumbens_basal_ganglia | rs843351 | C | 0.7139708399772644 | risk allele C increases intron usage of chr3:184276705-184276892(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **All switched-motif RBPs (recurrence across regions):** AGO2(2), CSTF2(2), ELAVL3(2), HNRNPD(2), HNRNPK(2), HNRNPU(2), KHDRBS3(2), PPIE(2), RBM24(2), RBM25(2), SRSF4(2), TRA2A(2), U2AF2(2)

## 5. Interpretation
ECE2 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
