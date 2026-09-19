# ANKRD36B — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair)

- **Traits:** SCZ  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.05 (Brain_Cerebellum)
- **Constraint:** LOEUF 0.569  ·  missense o/e 0.83

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| scz | sQTL | Brain_Cerebellum | 0.05 | no | — | no | — | direction unresolved |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Cerebellum | rs11695197 | G | nan | unresolved (variant not in GTEx signif_pairs or allele mismatch) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** G3BP2, HNRNPLL, PABPC3, IGF2BP2, YTHDC1, ADAR
- **All switched-motif RBPs (recurrence across regions):** ACO1(1), ADAR(1), AGO1(1), AGO2(1), CELF6(1), CNOT4(1), CPEB1(1), CPEB4(1), CSTF2(1), DDX19B(1), DHX58(1), EIF4A3(1), ESRP1(1), G3BP2(1), HNRNPA1L2(1), HNRNPAB(1), HNRNPCL1(1), HNRNPDL(1), HNRNPLL(1), HNRNPM(1)

## 5. Interpretation
ANKRD36B has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.
