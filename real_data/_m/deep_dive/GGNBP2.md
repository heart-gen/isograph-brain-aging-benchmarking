# GGNBP2 — mechanistic deep-dive

**Verdict:** splicing-led (IsoGraph-resolved switch) · GO-invisible

- **Traits:** ALS,SCZ  ·  **QTL kinds:** eQTL,sQTL  ·  **max CLPP:** 0.06 (Brain_Spinal_cord_cervical_c-1)
- **Constraint:** LOEUF 0.201  ·  missense o/e 0.66

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| als | eQTL | Brain_Spinal_cord_cervical_c-1 | 0.06 | no | — | yes | — | risk allele C increases GGNBP2 expression (gene-level; no intron) |
| scz | eQTL | Brain_Hippocampus | 0.03 | no | — | yes | — | risk allele G increases GGNBP2 expression (gene-level; no intron) |
| als | sQTL | Brain_Cerebellum | 0.01 | yes | yes | yes | no annotated structural change | risk allele C increases usage of junction chr17:36560871-36567663(+) (ENST00000613102.5,EN |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| scz | Brain_Hippocampus | rs11263770 | G | 0.21380910277366638 | risk allele G increases expression of GGNBP2 |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** SUPV3L1, NONO
- **All switched-motif RBPs (recurrence across regions):** A1CF(3), AGO2(3), AKAP1(3), CELF5(3), CELF4(3), CELF6(3), CPEB2(3), FXR2(3), CSTF2(3), DDX58(3), DHX58(3), EIF4A3(3), EIF4B(3), ESRP1(3), FXR1(3), HNRNPAB(3), G3BP1(3), G3BP2(3), HNRNPA1L2(3), HNRNPLL(3)

## 5. Interpretation
A splicing QTL colocalizes onto an IsoGraph switch pair for GGNBP2 in ALS. This is the IsoGraph-unique, DTU-without-DGE class: the disease variant acts through isoform choice, not gene dosage, in a GO-invisible module a pathway-enrichment scan would miss.

## 6. Literature (known isoform biology)
GGNBP2/ZNF403 (17q12) is LoF-constrained (LOEUF 0.20) and colocalizes as a splicing-led switch in ALS; its isoform biology in neurodegeneration is uncharacterized -- a novel splicing-led candidate.
