# MYO15A — mechanistic deep-dive

**Verdict:** expression-led (eQTL gene-level) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** eQTL  ·  **max CLPP:** 0.02 (Brain_Amygdala)
- **Constraint:** LOEUF 0.889  ·  missense o/e 0.99

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | eQTL | Brain_Amygdala | 0.02 | no | — | yes | — | risk allele T increases MYO15A expression (gene-level; no intron) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** YTHDC1, DHX58, HNRNPAB, AKAP1, RBMS3, PABPC4, ENOX1, EIF4B, ESRP1
- **All switched-motif RBPs (recurrence across regions):** AKAP1(4), CPEB1(4), CPEB4(4), DDX19B(4), DHX58(4), EIF4B(4), ELAVL3(4), ENOX1(4), ESRP1(4), G3BP1(4), HNRNPA1L2(4), HNRNPAB(4), KHDRBS1(4), KHDRBS2(4), PABPC4(4), PTBP2(4), TRA2A(4), RBM25(4), RBMS3(4), SAMD4A(4)

## 5. Interpretation
MYO15A colocalizes as an eQTL (gene-level expression), with no splicing event resolving to an IsoGraph switch pair — an honest expression-confounded case that abundance networks would also capture.
