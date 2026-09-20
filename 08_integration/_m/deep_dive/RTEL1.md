# RTEL1 — mechanistic deep-dive

**Verdict:** splicing (sQTL not resolved to switch pair) · GO-invisible

- **Traits:** AD  ·  **QTL kinds:** sQTL  ·  **max CLPP:** 0.01 (Brain_Cerebellar_Hemisphere)
- **Constraint:** LOEUF 0.638  ·  missense o/e 1.10

## 1–3. Genetic anchor → switch → coding consequence
| trait | kind | tissue | CLPP | in switch pair | concordant | GO-inv | struct. consequence | resolved event |
|-------|------|--------|------|----------------|------------|--------|---------------------|----------------|
| ad | sQTL | Brain_Cerebellar_Hemisphere | 0.01 | no | — | yes | — | risk allele C increases usage of junction chr20:63661496-63661762(+) (unmapped transcript; |

### Signed risk-allele direction (colocalized loci with allele matching)
| trait | tissue | rsID | risk allele | risk QTL effect | direction |
|-------|--------|------|-------------|-----------------|-----------|
| ad | Brain_Cerebellar_Hemisphere | rs61753459 | C | 0.5883097648620605 | risk allele C increases intron usage of chr20:63661496-63661762(+) |

## 4. Regulatory logic (RBP motifs in switched exons)
- **Switched *and* module-enriched (q<0.05) RBPs:** RBM6, IGF2BP1, HNRNPA3, CELF6, YTHDC1, ACO1, ZFP36L2, PABPC4, RBM24, RALY, SART3, RBM25, CSTF2, SRSF10, NELFE
- **All switched-motif RBPs (recurrence across regions):** DDX19B(5), CSTF2(5), CPEB4(5), ELAVL3(5), EIF4A3(5), ELAVL4(5), DHX58(5), HNRNPCL1(5), HNRNPC(5), HNRNPA0(5), ESRP1(5), YTHDC1(5), U2AF2(5), SYNCRIP(5), TRA2A(5), HNRNPU(5), HNRNPD(5), KHDRBS2(5), IGHMBP2(5), PABPC4(5)

## 5. Interpretation
RTEL1 has a colocalizing sQTL, but the junction does not map onto the IsoGraph switch pair for the tissue — splicing-associated but not resolved to a switch; a candidate for deeper transcript-level follow-up.

## 6. Literature (known isoform biology)
RTEL1 (telomere-maintenance helicase; AD/SCZ locus) has documented alternative C-terminal isoforms in other tissues, but a brain disease-specific splice role is not established. Not in the current anchored set.

_Curation: novel_candidate._ No established disease-specific isoform biology was found. That is a statement about the literature, not about the evidence here: an under-characterised switch is what this method is built to surface, and this row is a nomination rather than a null result.
